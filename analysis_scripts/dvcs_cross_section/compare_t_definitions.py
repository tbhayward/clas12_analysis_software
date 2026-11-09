#!/usr/bin/env python3
"""Quick DVCS candidate count migration: proton |t| versus photon |t|."""
import argparse
from array import array
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EDGES = (0.11, 0.15, 0.25, 0.40, 0.60, 0.80, 1.00)
LABELS = [f"{a:.2f}–{b:.2f}" for a,b in zip(EDGES[:-1], EDGES[1:])]
MIGRATION_EDGES = (-1.0e9, 0.0, *EDGES, 1.0e9)
MIGRATION_LABELS = ("unphysical", "<0.11", *LABELS, ">=1.00")
DELTA_EDGES = np.linspace(-1.0, 3.0, 201)
BEAM = {"fa18_inb":10.604,"fa18_out":10.604,
        "sp18_inb":10.594,"sp18_out":10.594,"sp19_inb":10.2}
MC_CURRENT = {"fa18_inb":"50nA","fa18_out":"50nA",
              "sp18_inb":"50nA","sp18_out":"45nA","sp19_inb":"50nA"}

def binning(merge_low_t=False):
    edges = (0.11, 0.25, 0.40, 0.60, 0.80, 1.00) if merge_low_t else EDGES
    labels = [f"{a:.2f}–{b:.2f}" for a,b in zip(edges[:-1],edges[1:])]
    migration_edges = (-1.0e9, 0.0, *edges, 1.0e9)
    migration_labels = ("unphysical", "<0.11", *labels, ">=1.00")
    return edges, labels, migration_edges, migration_labels

def output_path(root,period,merge_low_t=False):
    return root/period/("merged_low_t" if merge_low_t else "")

def exclusivity_expression(path, period, branches, sample):
    """Use the data or MC windows and topology keys read by total_counts.cpp."""
    if not path.is_file():
        raise FileNotFoundError(f"Production exclusivity cuts not found: {path}")
    document = json.loads(path.read_text())
    period_name = {"fa18":"Fa18","sp18":"Sp18","sp19":"Sp19"}[period[:4]]
    bend_name = "Inb" if period.endswith("inb") else "Out"
    topologies = ((1,1,"FD_FD"),(2,1,"CD_FD"),(2,0,"CD_FT"))
    allowed = {"Mx2","pTmiss","theta","theta_gamma_gamma","Delta_phi","Emiss2","Mx2_2"}
    expressions = []
    for detector1,detector2,topology in topologies:
        key = f"DVCS_{period_name}_{bend_name}_{topology}"
        if key not in document or sample not in document[key]:
            raise ValueError(f"Missing {sample} exclusivity cuts for {key}")
        block = document[key][sample]
        if "Mx2" not in block or not set(block).issubset(allowed):
            raise ValueError(f"Invalid production cut variables for {key}: {sorted(block)}")
        conditions = [f"detector1=={detector1}",f"detector2=={detector2}"]
        for variable,stats in block.items():
            if variable == "Delta_phi" and variable not in branches:
                if not {"p1_phi","p2_phi"}.issubset(branches):
                    raise ValueError("Delta_phi and fallback phi branches missing")
                value = "std::abs(std::remainder(p1_phi-p2_phi,6.283185307179586))"
            else:
                if variable not in branches:
                    raise ValueError(f"Missing cut branch {variable} for {key}")
                value = variable
            high,low = stats.get("cut_high"),stats.get("cut_low")
            if stats.get("mode") == "upper_quantile":
                if high is None or not math.isfinite(float(high)):
                    raise ValueError(f"Missing upper cut for {key}/{variable}")
                conditions += [f"std::isfinite({value})",f"{value}<={float(high):.17g}"]
            else:
                if (low is None or high is None or not math.isfinite(float(low))
                        or not math.isfinite(float(high)) or float(high)<=float(low)):
                    mean,std = float(stats["mean"]),float(stats["std"])
                    low,high = mean-3*std,mean+3*std
                conditions += [f"std::isfinite({value})",
                               f"{value}>={float(low):.17g}",
                               f"{value}<={float(high):.17g}"]
        expressions.append("("+" && ".join(conditions)+")")
        print(f"[EXCLUSIVITY] {sample} {key}: {', '.join(block)}")
    return " || ".join(expressions)

def plot_summary(data, mc, target, period):
    n = len(LABELS)
    nm = len(MIGRATION_LABELS)
    samples=[("Data",data)] + ([("Reconstructed MC",mc)] if mc is not None else [])
    fig,axes=plt.subplots(len(samples),3,figsize=(19,5.5*len(samples)),
                          squeeze=False,gridspec_kw={"width_ratios":[.85,1.4,1.1]},
                          constrained_layout=True)
    for row,(sample,values) in enumerate(samples):
        proton,photon,migration,delta_counts,candidates=values
        ax1,ax2,ax3=axes[row]
        ratio=np.divide(photon,proton,out=np.full(n,np.nan),where=proton>0)
        matrix=migration[2:2+n,:].astype(float)
        row_sum=matrix.sum(axis=1,keepdims=True)
        percent=np.divide(matrix,row_sum,out=np.zeros_like(matrix),where=row_sum>0)*100
        displayed=np.round(percent,1)
        for i in range(n):
            if row_sum[i,0]>0:
                largest=int(np.argmax(percent[i]))
                displayed[i,largest]+=np.round(100.0-displayed[i].sum(),1)
        xx=np.arange(n)
        ax1.bar(xx,ratio,color="#2673a8")
        ax1.axhline(1,color="#b42318",linestyle="--",linewidth=1.2)
        for x,value in zip(xx,ratio):
            if np.isfinite(value):
                ax1.text(x,value+.035,f"{value:.2f}",ha="center",fontsize=9)
        finite=ratio[np.isfinite(ratio)]
        ymax=max(1.55,float(finite.max())*1.15) if len(finite) else 1.55
        ax1.set(xticks=xx,xticklabels=LABELS,ylim=(0,ymax),
                xlabel="|t| bin (GeV²)",ylabel="Photon / proton count",
                title=f"{sample}: count ratio (N={candidates:,})")
        ax1.tick_params(axis="x",rotation=30)
        image=ax2.imshow(percent,origin="lower",aspect="auto",cmap="Blues",vmin=0,vmax=100)
        ax2.set(xticks=np.arange(nm),yticks=np.arange(n),
                xticklabels=MIGRATION_LABELS,yticklabels=LABELS,
                xlabel="Photon-defined |t| bin (GeV²)",
                ylabel="Proton-defined |t| bin (GeV²)",
                title=f"{sample}: destination (each row = 100.0%)")
        ax2.tick_params(axis="x",rotation=35)
        for i in range(n):
            for j in range(nm):
                if displayed[i,j]>=0.1:
                    ax2.text(j,i,f"{displayed[i,j]:.1f}%",ha="center",va="center",
                             fontsize=7.5,color="white" if percent[i,j]>=45 else "#162b45")
        fig.colorbar(image,ax=ax2,label="Percent of proton-defined bin")
        if delta_counts is not None:
            colors=plt.cm.viridis(np.linspace(.08,.9,6))
            for label,counts,color in zip(LABELS,delta_counts,colors):
                total=counts.sum()
                if total:
                    ax3.stairs(counts/total,DELTA_EDGES,label=label,color=color,linewidth=1.6)
            ax3.axvline(0,color="#b42318",linestyle="--",linewidth=1.1)
            ax3.legend(title="Proton |t| (GeV²)",fontsize=7,title_fontsize=8,
                       loc="upper right",frameon=True)
        else:
            ax3.text(.5,.5,"Run event scan to populate Δ|t|",ha="center",va="center",
                     transform=ax3.transAxes)
        ax3.set(xlim=(-.6,1.0),xlabel="Δ|t| = |t| photon − |t| proton (GeV²)",
                ylabel="Fraction per 0.02 GeV²",title=f"{sample}: event-level shift")
        ax3.grid(alpha=.15)
    fig.suptitle(f"{period}: proton and photon reconstructed |t| comparison")
    target.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(target,dpi=170)
    plt.close(fig)
    print(f"[PLOT] {target}")

def result_path(out,sample,name):
    return out/(name if sample=="data" else "mc_"+name)

def read_result(out,sample):
    counts_path=result_path(out,sample,"counts.csv")
    migration_path=result_path(out,sample,"migration.csv")
    if not counts_path.is_file() or not migration_path.is_file():
        return None
    with counts_path.open() as f:
        rows=list(csv.DictReader(f))
    proton=np.array([int(row["proton_count"]) for row in rows])
    photon=np.array([int(row["photon_count"]) for row in rows])
    with migration_path.open() as f:
        matrix=list(csv.reader(f))[1:]
    migration=np.array([[int(value) for value in row[1:]] for row in matrix])
    delta_counts=None
    delta_path=result_path(out,sample,"delta_t_hist.csv")
    if delta_path.is_file():
        with delta_path.open() as f:
            delta_rows=list(csv.DictReader(f))
        delta_counts=np.array([[int(row[f"proton_bin_{i+1}"])
                                for row in delta_rows] for i in range(len(LABELS))])
    return proton,photon,migration,delta_counts,int(migration.sum())

def plot_available(out,period):
    data=read_result(out,"data")
    mc=read_result(out,"mc")
    if data is None:
        print(f"[PLOT] waiting for data CSVs in {out}")
        return
    plot_summary(data,mc,out/"t_definition_counts.png",period)
    if mc is not None:
        d=np.divide(data[1],data[0],out=np.full(len(LABELS),np.nan),where=data[0]>0)
        r=np.divide(mc[1],mc[0],out=np.full(len(LABELS),np.nan),where=mc[0]>0)
        print("[DIAGNOSTIC] data-ratio / MC-ratio by |t| bin (candidate counts, not cross sections)")
        for label,dr,mr in zip(LABELS,d,r):
            print(f"  {label}: data={dr:.3f} MC={mr:.3f} double_ratio={dr/mr:.3f}")

def main():
    global EDGES, LABELS, MIGRATION_EDGES, MIGRATION_LABELS
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-root", type=Path,
        default=Path("/volatile/clas12/thayward/dvcs_pass2_full_rebuild"))
    ap.add_argument("--period", choices=BEAM, default="fa18_inb")
    ap.add_argument("--max-files", type=int, default=1,
                    help="Use 0 for all parts")
    ap.add_argument("--skip-files", type=int, default=0,
                    help="Skip this many sorted parts before selecting files")
    ap.add_argument("--output-dir", type=Path,
                    default=Path("output/t_definition_quick_test"))
    ap.add_argument("--extra-cut", default="1",
                    help="Optional additional ROOT selection expression")
    ap.add_argument("--cuts-json", type=Path,
                    default=Path(__file__).resolve().parent/"output/jsons/combined_cuts.json",
                    help="Saved production exclusivity windows")
    ap.add_argument("--plot-only",action="store_true",
                    help="Replot existing counts.csv and migration.csv without reading ROOT")
    ap.add_argument("--sample",choices=("data","mc"),default="data",
                    help="Sample to scan; plot includes both rows once both CSVs exist")
    ap.add_argument("--workers",type=int,default=10,
                    help="ROOT processing threads; maximum 10")
    ap.add_argument("--merge-low-t",action="store_true",
                    help="Combine 0.11–0.15 and 0.15–0.25 into 0.11–0.25 for this run")
    args = ap.parse_args()
    EDGES, LABELS, MIGRATION_EDGES, MIGRATION_LABELS = binning(args.merge_low_t)
    out = output_path(args.output_dir,args.period,args.merge_low_t)
    if args.plot_only:
        plot_available(out,args.period)
        return
    if not 1<=args.workers<=10:
        ap.error("--workers must be between 1 and 10")
    import ROOT
    ROOT.EnableImplicitMT(args.workers)
    if args.sample=="data":
        parts=args.input_root/"data"/args.period/"dvcs_data"/"parts"
        paths=sorted(parts.glob("*.root"))
    else:
        parts=args.input_root/"mc"/args.period/"dvcs_mc"/MC_CURRENT[args.period]/"parts"
        paths=sorted(parts.glob("*.rec_only.root"))
    if not paths:
        ap.error(f"No {args.sample} PhysicsEvents parts found in {parts}")
    paths = paths[args.skip_files:]
    if args.max_files > 0:
        paths = paths[:args.max_files]
    if not paths:
        ap.error("No parts left after --skip-files")
    chain = ROOT.TChain("PhysicsEvents")
    for path in paths:
        chain.Add(str(path))
    needed = {"t1","e_p","e_theta","e_phi","p2_p","p2_theta","p2_phi",
              "Q2","x","detector1","detector2"}
    branches = {b.GetName() for b in chain.GetListOfBranches()}
    missing = needed - branches
    if missing:
        ap.error(f"Missing branches: {sorted(missing)}")
    try:
        exclusive = exclusivity_expression(args.cuts_json,args.period,branches,args.sample)
    except (FileNotFoundError, ValueError) as exc:
        ap.error(str(exc))
    ROOT.gInterpreter.Declare(r'''
    #include <cmath>
    double quick_tgamma(double ep,double eth,double eph,double gp,double gth,double gph,double E) {
      constexpr double me=0.00051099895;
      double ex=ep*std::sin(eth)*std::cos(eph), ey=ep*std::sin(eth)*std::sin(eph), ez=ep*std::cos(eth);
      double gx=gp*std::sin(gth)*std::cos(gph), gy=gp*std::sin(gth)*std::sin(gph), gz=gp*std::cos(gth);
      double q0=E-std::sqrt(ep*ep+me*me)-gp;
      double qx=-ex-gx, qy=-ey-gy, qz=E-ez-gz;
      return q0*q0-qx*qx-qy*qy-qz*qz;
    }''')
    # Hold the proton-defined global |t| range fixed, then switch only bin assignment.
    stages = (
        ("all entries", "1"),
        ("physical proton t", "std::isfinite(t1) && t1<0"),
        ("proton |t| < 1", "-t1<1.0"),
        ("Q2 > 1", "std::isfinite(Q2) && Q2>1"),
        ("valid x", "std::isfinite(x) && x>0 && x<1"),
        ("W > 2", "0.9382720813*0.9382720813+Q2*(1/x-1)>4"),
        ("electron p > 2", "std::isfinite(e_p) && e_p>2"),
        ("photon p > 0", "std::isfinite(p2_p) && p2_p>0"),
        ("DVCS topology", "(detector1==1 || detector1==2) && (detector2==0 || detector2==1)"),
        ("saved exclusivity cuts", exclusive),
        ("extra cut", args.extra_cut),
    )
    df = ROOT.RDataFrame(chain)
    cutflow = []
    for name, condition in stages:
        df = df.Filter(condition)
        cutflow.append((name, df.Count()))
    df = df.Define("tp", "-t1").Define(
        "t_gamma_signed",f"quick_tgamma(e_p,e_theta,e_phi,p2_p,p2_theta,p2_phi,{BEAM[args.period]:.17g})")
    df = df.Filter("std::isfinite(t_gamma_signed)")
    nonphysical = df.Filter("t_gamma_signed>=0").Count()
    df = df.Define("tg", "-t_gamma_signed")
    df = df.Define("delta_t", "tg-tp")
    edges = array("d",EDGES)
    n = len(LABELS)
    nm = len(MIGRATION_LABELS)
    hp = df.Histo1D(ROOT.RDF.TH1DModel("hp","",n,edges),"tp")
    hg = df.Histo1D(ROOT.RDF.TH1DModel("hg","",n,edges),"tg")
    medges = array("d",MIGRATION_EDGES)
    hm = df.Histo2D(ROOT.RDF.TH2DModel("hm","",nm,medges,nm,medges),"tp","tg")
    delta_handles=[]
    delta_edges=array("d",DELTA_EDGES)
    for i,(low,high) in enumerate(zip(EDGES[:-1],EDGES[1:])):
        cell=df.Filter(f"tp>={low:.17g} && tp<{high:.17g} && tg>0")
        delta_handles.append(cell.Histo1D(
            ROOT.RDF.TH1DModel(f"hdelta_{i}","",len(DELTA_EDGES)-1,delta_edges),
            "delta_t"))
    count = df.Count()
    ROOT.RDF.RunGraphs([hp,hg,hm,count,nonphysical,*delta_handles,
                        *[handle for _,handle in cutflow]])
    print("[INPUT]", *[str(p) for p in paths], sep="\n  ")
    for name, handle in cutflow:
        print(f"[CUTFLOW] {name}: {int(handle.GetValue()):,}")
    a = np.array([hp.GetValue().GetBinContent(i+1) for i in range(n)],dtype=int)
    b = np.array([hg.GetValue().GetBinContent(i+1) for i in range(n)],dtype=int)
    migration = np.array([[hm.GetValue().GetBinContent(i+1,j+1)
                           for j in range(nm)] for i in range(nm)],dtype=int)
    delta_counts=np.array([[handle.GetValue().GetBinContent(j+1)
                            for j in range(len(DELTA_EDGES)-1)]
                           for handle in delta_handles],dtype=int)
    out.mkdir(parents=True,exist_ok=True)
    with result_path(out,args.sample,"counts.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["t_min_GeV2","t_max_GeV2","proton_count","photon_count",
                    "photon_over_proton"])
        for i in range(n):
            w.writerow([EDGES[i],EDGES[i+1],a[i],b[i],b[i]/a[i] if a[i] else ""])
    with result_path(out,args.sample,"migration.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["proton_bin / photon_bin",*MIGRATION_LABELS])
        for label,row in zip(MIGRATION_LABELS,migration):
            w.writerow([label,*row])
    with result_path(out,args.sample,"delta_t_hist.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["delta_t_min_GeV2","delta_t_max_GeV2",
                    *[f"proton_bin_{i+1}" for i in range(n)]])
        for j in range(len(DELTA_EDGES)-1):
            w.writerow([DELTA_EDGES[j],DELTA_EDGES[j+1],
                        *delta_counts[:,j]])
    plot_available(out,args.period)
    print(f"[SELECTION] baseline Q²>1, W>2, electron p>2, proton |t|<1; saved {args.sample} exclusivity cuts")
    print(f"[SAMPLE] {args.sample} [FILES] {len(paths)} [CANDIDATES] {int(count.GetValue()):,}")
    print(f"[PHOTON T] unphysical signed t_gamma >= 0: {int(nonphysical.GetValue()):,}")
    for label, n1, n2 in zip(LABELS,a,b):
        print(f"|t| {label}: proton={n1:,} photon={n2:,} ratio={n2/n1:.3f}" if n1 else
              f"|t| {label}: proton=0 photon={n2:,}")
    print(f"[LOWEST BIN] stays={migration[2,2]:,} leaves_to_unphysical={migration[2,0]:,} "
          f"leaves_to_below={migration[2,1]:,} leaves_to_higher={migration[2,3:].sum():,} "
          f"enters_from_below={migration[1,2]:,} enters_from_higher={migration[3:,2].sum():,}")
    print(f"[OUTPUT] {out}")

if __name__=="__main__":
    main()
