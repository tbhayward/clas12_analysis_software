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
BEAM = {"fa18_inb":10.604,"fa18_out":10.604,
        "sp18_inb":10.594,"sp18_out":10.594,"sp19_inb":10.2}

def exclusivity_expression(path, period, branches):
    """Use the data windows and topology keys read by total_counts.cpp."""
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
        if key not in document or "data" not in document[key]:
            raise ValueError(f"Missing data exclusivity cuts for {key}")
        block = document[key]["data"]
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
        print(f"[EXCLUSIVITY] {key}: {', '.join(block)}")
    return " || ".join(expressions)

def main():
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
                    default=Path("output/jsons/combined_cuts.json"),
                    help="Saved production exclusivity windows")
    args = ap.parse_args()
    import ROOT
    paths = sorted((args.input_root/"data"/args.period/"dvcs_data"/"parts").glob("*.root"))
    if not paths:
        ap.error("No DVCS data parts found")
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
        exclusive = exclusivity_expression(args.cuts_json,args.period,branches)
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
    edges = array("d",EDGES)
    hp = df.Histo1D(ROOT.RDF.TH1DModel("hp","",6,edges),"tp")
    hg = df.Histo1D(ROOT.RDF.TH1DModel("hg","",6,edges),"tg")
    medges = array("d",MIGRATION_EDGES)
    hm = df.Histo2D(ROOT.RDF.TH2DModel("hm","",9,medges,9,medges),"tp","tg")
    count = df.Count()
    ROOT.RDF.RunGraphs([hp,hg,hm,count,nonphysical,*[handle for _,handle in cutflow]])
    print("[INPUT]", *[str(p) for p in paths], sep="\n  ")
    for name, handle in cutflow:
        print(f"[CUTFLOW] {name}: {int(handle.GetValue()):,}")
    a = np.array([hp.GetValue().GetBinContent(i+1) for i in range(6)],dtype=int)
    b = np.array([hg.GetValue().GetBinContent(i+1) for i in range(6)],dtype=int)
    migration = np.array([[hm.GetValue().GetBinContent(i+1,j+1)
                           for j in range(9)] for i in range(9)],dtype=int)
    out = args.output_dir/args.period
    out.mkdir(parents=True,exist_ok=True)
    with (out/"counts.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["t_min_GeV2","t_max_GeV2","proton_count","photon_count",
                    "photon_over_proton"])
        for i in range(6):
            w.writerow([EDGES[i],EDGES[i+1],a[i],b[i],b[i]/a[i] if a[i] else ""])
    with (out/"migration.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["proton_bin / photon_bin",*MIGRATION_LABELS])
        for label,row in zip(MIGRATION_LABELS,migration):
            w.writerow([label,*row])
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4.7),constrained_layout=True)
    xx=np.arange(6)
    ax1.bar(xx-.19,a,width=.38,label="Proton definition")
    ax1.bar(xx+.19,b,width=.38,label="Photon definition")
    ax1.set(xticks=xx,xticklabels=LABELS,xlabel="|t| bin (GeV²)",
            ylabel="Candidate count",title="Same events; different |t| assignment")
    ax1.tick_params(axis="x",rotation=30)
    ax1.legend()
    image=ax2.imshow(migration,origin="lower",aspect="auto",cmap="Blues")
    mx=np.arange(9)
    ax2.set(xticks=mx,yticks=mx,xticklabels=MIGRATION_LABELS,yticklabels=MIGRATION_LABELS,
            xlabel="Photon |t| bin (GeV²)",ylabel="Proton |t| bin (GeV²)",
            title="Bin migration")
    ax2.tick_params(axis="x",rotation=30)
    fig.colorbar(image,ax=ax2,label="Candidates")
    fig.suptitle(f"{args.period}: {len(paths)} file(s), {int(count.GetValue()):,} baseline candidates")
    fig.savefig(out/"t_definition_counts.png",dpi=160)
    print(f"[SELECTION] baseline Q²>1, W>2, electron p>2, proton |t|<1; saved data exclusivity cuts")
    print(f"[FILES] {len(paths)} [CANDIDATES] {int(count.GetValue()):,}")
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
