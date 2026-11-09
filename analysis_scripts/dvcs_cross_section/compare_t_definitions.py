#!/usr/bin/env python3
"""Quick DVCS candidate count migration: proton |t| versus photon |t|."""
import argparse
from array import array
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EDGES = (0.11, 0.15, 0.25, 0.40, 0.60, 0.80, 1.00)
LABELS = [f"{a:.2f}–{b:.2f}" for a,b in zip(EDGES[:-1], EDGES[1:])]
MIGRATION_EDGES = (0.0, *EDGES, 2.0)
MIGRATION_LABELS = ("<0.11", *LABELS, "1.00–2.00")
BEAM = {"fa18_inb":10.604,"fa18_out":10.604,
        "sp18_inb":10.594,"sp18_out":10.594,"sp19_inb":10.2}

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-root", type=Path,
        default=Path("/volatile/clas12/thayward/dvcs_pass2_full_rebuild"))
    ap.add_argument("--period", choices=BEAM, default="fa18_inb")
    ap.add_argument("--max-files", type=int, default=1,
                    help="Use 0 for all parts")
    ap.add_argument("--output-dir", type=Path,
                    default=Path("output/t_definition_quick_test"))
    ap.add_argument("--extra-cut", default="1",
                    help="Optional additional ROOT selection expression")
    args = ap.parse_args()
    import ROOT
    paths = sorted((args.input_root/"data"/args.period/"dvcs_data"/"parts").glob("*.root"))
    if not paths:
        ap.error("No DVCS data parts found")
    if args.max_files > 0:
        paths = paths[:args.max_files]
    chain = ROOT.TChain("PhysicsEvents")
    for path in paths:
        chain.Add(str(path))
    needed = {"t1","e_p","e_theta","e_phi","p2_p","p2_theta","p2_phi",
              "Q2","x","detector1","detector2"}
    missing = needed - {b.GetName() for b in chain.GetListOfBranches()}
    if missing:
        ap.error(f"Missing branches: {sorted(missing)}")
    ROOT.gInterpreter.Declare(r'''
    #include <cmath>
    double quick_tgamma(double ep,double eth,double eph,double gp,double gth,double gph,double E) {
      constexpr double me=0.00051099895;
      double ex=ep*std::sin(eth)*std::cos(eph), ey=ep*std::sin(eth)*std::sin(eph), ez=ep*std::cos(eth);
      double gx=gp*std::sin(gth)*std::cos(gph), gy=gp*std::sin(gth)*std::sin(gph), gz=gp*std::cos(gth);
      double q0=E-std::sqrt(ep*ep+me*me)-gp;
      double qx=-ex-gx, qy=-ey-gy, qz=E-ez-gz;
      return std::abs(q0*q0-qx*qx-qy*qy-qz);
    }''')
    # Hold the proton-defined global |t| range fixed, then switch only bin assignment.
    base = ("std::isfinite(t1) && t1<0 && -t1<1.0 && "
            "std::isfinite(Q2) && Q2>1 && std::isfinite(x) && x>0 && x<1 && "
            "0.9382720813*0.9382720813+Q2*(1/x-1)>4 && "
            "std::isfinite(e_p) && e_p>2 && "
            "std::isfinite(p2_p) && p2_p>0 && "
            "(detector1==1 || detector1==2) && "
            "(detector2==0 || detector2==1)")
    df = ROOT.RDataFrame(chain).Filter(base).Filter(args.extra_cut)
    df = df.Define("tp", "-t1").Define(
        "tg",f"quick_tgamma(e_p,e_theta,e_phi,p2_p,p2_theta,p2_phi,{BEAM[args.period]:.17g})")
    df = df.Filter("std::isfinite(tg)")
    edges = array("d",EDGES)
    hp = df.Histo1D(ROOT.RDF.TH1DModel("hp","",6,edges),"tp")
    hg = df.Histo1D(ROOT.RDF.TH1DModel("hg","",6,edges),"tg")
    medges = array("d",MIGRATION_EDGES)
    hm = df.Histo2D(ROOT.RDF.TH2DModel("hm","",8,medges,8,medges),"tp","tg")
    count = df.Count()
    ROOT.RDF.RunGraphs([hp,hg,hm,count])
    a = np.array([hp.GetValue().GetBinContent(i+1) for i in range(6)],dtype=int)
    b = np.array([hg.GetValue().GetBinContent(i+1) for i in range(6)],dtype=int)
    migration = np.array([[hm.GetValue().GetBinContent(i+1,j+1)
                           for j in range(8)] for i in range(8)],dtype=int)
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
    mx=np.arange(8)
    ax2.set(xticks=mx,yticks=mx,xticklabels=MIGRATION_LABELS,yticklabels=MIGRATION_LABELS,
            xlabel="Photon |t| bin (GeV²)",ylabel="Proton |t| bin (GeV²)",
            title="Bin migration")
    ax2.tick_params(axis="x",rotation=30)
    fig.colorbar(image,ax=ax2,label="Candidates")
    fig.suptitle(f"{args.period}: {len(paths)} file(s), {int(count.GetValue()):,} baseline candidates")
    fig.savefig(out/"t_definition_counts.png",dpi=160)
    print(f"[SELECTION] baseline Q²>1, W>2, electron p>2, proton |t|<1; no exclusivity cuts")
    print(f"[FILES] {len(paths)} [CANDIDATES] {int(count.GetValue()):,}")
    for label, n1, n2 in zip(LABELS,a,b):
        print(f"|t| {label}: proton={n1:,} photon={n2:,} ratio={n2/n1:.3f}" if n1 else
              f"|t| {label}: proton=0 photon={n2:,}")
    print(f"[LOWEST BIN] stays={migration[1,1]:,} leaves_to_below={migration[1,0]:,} "
          f"leaves_to_higher={migration[1,2:].sum():,} "
          f"enters_from_below={migration[0,1]:,} "
          f"enters_from_higher={migration[2:,1].sum():,}")
    print(f"[OUTPUT] {out}")

if __name__=="__main__":
    main()
