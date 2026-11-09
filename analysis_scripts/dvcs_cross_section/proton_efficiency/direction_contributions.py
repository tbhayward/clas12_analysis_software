#!/usr/bin/env python3
"""Diagnose which charged track shifts the predicted missing-proton angle.

This reads completed MC Tags trees from the proton-efficiency campaign. It uses
the same truth-near selection as kinematic_migration.py. For each event, it
replaces the truth electron and pion three-momenta by their reconstructed
values in all possible orders. Shapley averaging assigns an exact, additive
bookkeeping share of the *polar-angle* shift to each charged track. The
allocation convention does not by itself identify a detector calibration
problem.
"""

import argparse
from array import array
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import math
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


MASS = .9382720813
THETA = np.array([0, 10, 20, 30, 40, 50, 60, 70], dtype=float)
P = np.array([.2, .5, .7, .9, math.sqrt(1 + 1/(4*MASS*MASS))])
TRACKS = ("e", "pip", "pim")
FIELDS = ("e_shift", "pip_shift", "pim_shift", "total_shift")
REQUIRED = {
    "truth_tag_match_valid", "truth_missing_M2", "truth_missing_theta",
    "truth_missing_p", "truth_missing_px", "truth_missing_py",
    "truth_missing_pz", "truth_best_index", "truth_best_match_angle",
    "truth_best_p", "truth_best_px", "truth_best_py", "truth_best_pz",
    "missing_theta", "missing_p", "missing_M2", "missing_px", "missing_py",
    "missing_pz", "e_p", "truth_e_p", "e_theta", "beam_energy", "Q2",
}
for track in TRACKS:
    REQUIRED.update(f"{track}_{part}" for part in ("p", "theta", "phi"))
    REQUIRED.update(f"truth_{track}_{part}" for part in ("p", "theta", "phi"))

CPP = r"""
#include <ROOT/RVec.hxx>
#include <array>
#include <cmath>
namespace proton_direction {
struct V3 { double x,y,z; };
V3 operator-(V3 a,V3 b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
V3 particle(double p,double theta,double phi) {
  double s=std::sin(theta);
  return {p*s*std::cos(phi),p*s*std::sin(phi),p*std::cos(theta)};
}
double polar(V3 v) {
  return std::atan2(std::hypot(v.x,v.y),v.z)*180.0/3.14159265358979323846;
}
ROOT::VecOps::RVec<double> split(
  double ep,double et,double ef,double tep,double tet,double tef,
  double pp,double pt,double pf,double tpp,double tpt,double tpf,
  double mp,double mt,double mf,double tmp,double tmt,double tmf,
  double tx,double ty,double tz) {
  constexpr double rad=3.14159265358979323846/180.0;
  // Reconstructed track angles are radians; truth track angles are degrees.
  std::array<V3,3> d{{
    particle(ep,et,ef)-particle(tep,tet*rad,tef*rad),
    particle(pp,pt,pf)-particle(tpp,tpt*rad,tpf*rad),
    particle(mp,mt,mf)-particle(tmp,tmt*rad,tmf*rad)}};
  V3 truth{tx,ty,tz};
  std::array<double,8> angle{};
  for (int mask=0;mask<8;++mask) {
    V3 v=truth;
    for (int i=0;i<3;++i) if (mask&(1<<i)) v=v-d[i];
    angle[mask]=polar(v);
  }
  ROOT::VecOps::RVec<double> result(4);
  for (int i=0;i<3;++i) {
    double value=0;
    for (int mask=0;mask<8;++mask) if (!(mask&(1<<i))) {
      int size=((mask&1)!=0)+((mask&2)!=0)+((mask&4)!=0);
      double weight=size==1 ? 1.0/6.0 : 1.0/3.0;
      value+=weight*(angle[mask|(1<<i)]-angle[mask]);
    }
    result[i]=value;
  }
  result[3]=angle[7]-angle[0];
  return result;
}
}
"""
SPLIT_ARGS = [value for track in TRACKS
              for value in (f"{track}_p", f"{track}_theta", f"{track}_phi",
                            f"truth_{track}_p", f"truth_{track}_theta",
                            f"truth_{track}_phi")]
SPLIT_ARGS += ["truth_missing_px", "truth_missing_py", "truth_missing_pz"]
_CPP_DECLARED = False


def completed_tasks(input_root, period, max_tags, match_angle):
    tasks, remaining = [], max_tags
    for path in sorted((input_root/period/"mc"/"parts").glob("*.root")):
        marker = path.with_suffix(".ok.json")
        if not marker.is_file():
            continue
        stop = None
        if max_tags:
            n = int(json.loads(marker.read_text()).get("selected_tags", 0))
            if n <= 0:
                continue
            stop = min(n, remaining)
            remaining -= stop
        tasks.append((str(path), stop, match_angle))
        if max_tags and remaining == 0:
            break
    if not tasks or (max_tags and remaining):
        raise RuntimeError("No complete MC parts or fewer tags than requested")
    return tasks


def read_part(task):
    global _CPP_DECLARED
    filename, stop, match_angle = task
    import ROOT

    ROOT.TH1.AddDirectory(False)
    source = ROOT.TFile.Open(filename, "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")
    try:
        tree = source.Get("Tags")
        if not tree:
            raise RuntimeError(f"Missing Tags tree in {filename}")
        branches = {b.GetName() for b in tree.GetListOfBranches()}
        if REQUIRED - branches:
            raise RuntimeError(f"Missing branches: {sorted(REQUIRED - branches)}")
        if not _CPP_DECLARED:
            ROOT.gInterpreter.Declare(CPP)
            _CPP_DECLARED = True
        frame = ROOT.RDataFrame(tree)
        if stop is not None:
            frame = frame.Range(0, min(stop, int(tree.GetEntries())))
        frame = frame.Filter(
            "truth_tag_match_valid > 0.5 && std::isfinite(truth_missing_M2) && "
            f"truth_missing_M2 >= {MASS*MASS-.15:.17g} && "
            f"truth_missing_M2 < {MASS*MASS+.15:.17g} && "
            "truth_best_index >= 0 && std::isfinite(truth_best_match_angle) && "
            f"truth_best_match_angle < {match_angle:.17g} && "
            "std::isfinite(truth_best_px) && std::isfinite(truth_best_py) && "
            "std::isfinite(truth_best_pz) && std::isfinite(truth_best_p) && "
            "truth_best_p > 0 && "
            "std::isfinite(truth_missing_p) && truth_missing_p > 0 && "
            "std::isfinite(missing_p) && missing_p > 0 && "
            "std::isfinite(missing_M2) && std::isfinite(e_p) && "
            "std::isfinite(truth_e_p) && std::isfinite(e_theta) && "
            "std::isfinite(beam_energy) && std::isfinite(Q2) && "
            "std::isfinite(truth_missing_px) && std::isfinite(truth_missing_py) && "
            "std::isfinite(truth_missing_pz) && std::isfinite(missing_px) && "
            "std::isfinite(missing_py) && std::isfinite(missing_pz) && "
            "std::isfinite(truth_missing_theta) && truth_missing_theta >= 0 && "
            "truth_missing_theta < 180 && std::isfinite(missing_theta) && "
            f"missing_theta >= {THETA[0]:.17g} && missing_theta < {THETA[-1]:.17g} && "
            f"missing_p >= {P[0]:.17g} && missing_p < {P[-1]:.17g}")
        for track in TRACKS:
            frame = frame.Filter(" && ".join(
                f"std::isfinite({prefix}{track}_{part})"
                for prefix in ("", "truth_") for part in ("p", "theta", "phi")))
        frame = frame.Define(
            "generated_proton_theta",
            "std::atan2(std::hypot(truth_best_px,truth_best_py),"
            "truth_best_pz)*57.29577951308232")
        frame = frame.Define(
            "angle_share", "proton_direction::split(" + ",".join(SPLIT_ARGS) + ")")
        for i, field in enumerate(FIELDS):
            frame = frame.Define(field, f"angle_share[{i}]")
            frame = frame.Define("abs_" + field, f"std::abs({field})")
        frame = frame.Define(
            "closure", "std::abs(e_shift+pip_shift+pim_shift-"
            "(missing_theta-truth_missing_theta))")
        count = frame.Count()
        failures = frame.Filter("closure > 0.05 || !std::isfinite(closure)").Count()
        tbins, pbins = array("d", THETA), array("d", P)
        def model(name):
            return ROOT.RDF.TH2DModel(
                name, "", len(THETA)-1, tbins, len(P)-1, pbins)
        hcount = frame.Histo2D(model("h_count"), "missing_theta", "missing_p")
        htruth = frame.Histo3D(
            ROOT.RDF.TH3DModel("h_truth_proton_angle", "",
                               180, array("d", np.arange(181)),
                               len(THETA)-1, tbins, len(P)-1, pbins),
            "generated_proton_theta", "missing_theta", "missing_p")
        sums = {field: frame.Histo2D(model("h_" + field),
                                    "missing_theta", "missing_p", field)
                for field in FIELDS}
        absolute = {field: frame.Histo2D(model("h_abs_" + field),
                                        "missing_theta", "missing_p",
                                        "abs_" + field)
                    for field in FIELDS}
        ROOT.RDF.RunGraphs([count, failures, hcount, htruth,
                            *sums.values(), *absolute.values()])
        def matrix(handle):
            hist = handle.GetValue()
            return np.array(
                [[hist.GetBinContent(i+1,j+1) for j in range(len(P)-1)]
                 for i in range(len(THETA)-1)], dtype=float)
        truth_histogram = htruth.GetValue()
        return {"count": int(count.GetValue()), "closure_failed": int(failures.GetValue()),
                "n": matrix(hcount),
                "truth_proton_angle": np.array(
                    [[[truth_histogram.GetBinContent(i+1,j+1,k+1)
                       for k in range(len(P)-1)]
                      for j in range(len(THETA)-1)]
                     for i in range(180)], dtype=np.int64),
                "sum": {key: matrix(value) for key,value in sums.items()},
                "abs": {key: matrix(value) for key,value in absolute.items()}}
    finally:
        source.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--period", default="fa18_inb")
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--max-tags", type=int, default=0)
    parser.add_argument("--match-angle", type=float, default=10.)
    parser.add_argument("--output-root", type=Path, default=Path("output/proton_efficiency"))
    args = parser.parse_args()
    if not 1 <= args.workers <= 10 or args.max_tags < 0 or args.match_angle <= 0:
        parser.error("workers must be 1..10; max-tags nonnegative; match-angle positive")
    tasks = completed_tasks(args.input_root, args.period, args.max_tags, args.match_angle)
    output = (args.output_root/args.period/
              (f"direction_contributions_{args.max_tags}" if args.max_tags
               else "direction_contributions_full"))
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    totals = {"count": 0, "closure_failed": 0,
              "n": np.zeros((len(THETA)-1,len(P)-1), dtype=float),
              "truth_proton_angle": np.zeros(
                  (180,len(THETA)-1,len(P)-1), dtype=np.int64),
              "sum": {key: np.zeros((len(THETA)-1,len(P)-1), dtype=float)
                      for key in FIELDS},
              "abs": {key: np.zeros((len(THETA)-1,len(P)-1), dtype=float)
                      for key in FIELDS}}
    print(f"[MC] files={len(tasks)} workers={min(args.workers,len(tasks))}", flush=True)
    with ProcessPoolExecutor(max_workers=min(args.workers,len(tasks))) as pool:
        futures = {pool.submit(read_part, task): task[0] for task in tasks}
        for done, future in enumerate(as_completed(futures), 1):
            try:
                part = future.result()
            except Exception as error:
                raise RuntimeError(f"Failed {futures[future]}: {error}") from error
            totals["count"] += part["count"]
            totals["closure_failed"] += part["closure_failed"]
            totals["n"] += part["n"]
            totals["truth_proton_angle"] += part["truth_proton_angle"]
            for group in ("sum", "abs"):
                for key in FIELDS:
                    totals[group][key] += part[group][key]
            if done % 10 == 0 or done == len(tasks):
                print(f"[MC] {done}/{len(tasks)} ({100*done/len(tasks):.1f}%) "
                      f"truth-near={totals['count']:,}", flush=True)
    if totals["closure_failed"]:
        raise RuntimeError(f"[VALIDATION_FAILED] polar-angle decomposition closure: "
                           f"{totals['closure_failed']:,} tags")
    if int(round(totals["n"].sum())) != totals["count"]:
        raise RuntimeError("[VALIDATION_FAILED] theta-momentum cell counts do not close")
    if int(totals["truth_proton_angle"].sum()) != totals["count"]:
        raise RuntimeError("[VALIDATION_FAILED] generated-proton-angle counts do not close")
    split_sum = sum((totals["sum"][key] for key in FIELDS[:3]),
                    np.zeros_like(totals["n"]))
    if not np.allclose(split_sum, totals["sum"]["total_shift"],
                       rtol=1e-8, atol=1e-5):
        raise RuntimeError("[VALIDATION_FAILED] cell-level direction shares do not close")
    with (output/"1_direction_cells.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("theta_low_deg","theta_high_deg","p_low_GeV","p_high_GeV",
                         "truth_near_tags","generated_proton_theta_lt_50",
                         "generated_proton_theta_lt_60",
                         *(f"mean_{key}_deg" for key in FIELDS),
                         *(f"mean_abs_{key}_deg" for key in FIELDS)))
        for i in range(len(THETA)-1):
            for j in range(len(P)-1):
                n = totals["n"][i,j]
                writer.writerow((THETA[i],THETA[i+1],P[j],P[j+1],int(n),
                                 int(totals["truth_proton_angle"][:50,i,j].sum()),
                                 int(totals["truth_proton_angle"][:60,i,j].sum()),
                                 *(totals["sum"][key][i,j]/n if n else math.nan
                                   for key in FIELDS),
                                 *(totals["abs"][key][i,j]/n if n else math.nan
                                   for key in FIELDS)))
    bands = (4,5,6)
    fig, axes = plt.subplots(1, 3, figsize=(14,4.5), sharey=True,
                             constrained_layout=True)
    labels = ("electron", "π+", "π−")
    width = .22
    for ax,i in zip(axes,bands):
        n = totals["n"][i,:]
        x = np.arange(len(n))
        for index,key in enumerate(FIELDS[:3]):
            value = np.divide(totals["sum"][key][i,:],n,
                              out=np.full(len(n),np.nan),where=n>0)
            ax.bar(x+(index-1)*width,value,width,label=labels[index])
        ax.axhline(0,color="black",linewidth=.8)
        ax.set_xticks(x,[f"{P[j]:.2f}–{P[j+1]:.2f}" for j in range(len(n))],
                      rotation=35,ha="right")
        ax.set_title(f"Predicted θ {THETA[i]:.0f}–{THETA[i+1]:.0f}°")
        ax.set_xlabel("Predicted proton p [GeV]")
        ax.grid(axis="y",alpha=.25)
    axes[0].set_ylabel("Mean contribution to θ_rec − θ_truth [degrees]")
    axes[-1].legend(frameon=False)
    fig.suptitle("Fa18 Inb: exact per-track decomposition of missing-proton polar-angle shift")
    fig.savefig(output/"2_direction_contributions.png",dpi=170)
    fig,ax = plt.subplots(figsize=(8.5,6),constrained_layout=True)
    h = totals["truth_proton_angle"][:,4:7,:].sum(axis=2)
    image = ax.pcolormesh(THETA[4:8],np.arange(181),
                          np.ma.masked_less_equal(h,0),
                          norm=LogNorm(vmin=1,vmax=max(10,int(h.max()))),
                          cmap="viridis",shading="auto")
    ax.plot([40,70],[40,70],"w--",linewidth=1.3)
    ax.set(xlabel="Reconstructed predicted proton θ [degrees]",
           ylabel="Generated proton θ [degrees]",xlim=(40,70),ylim=(0,90),
           title="Generated proton angle versus reconstructed prediction")
    fig.colorbar(image,ax=ax,label="Truth-near tags / 1° × 10° cell")
    fig.savefig(output/"3_generated_proton_migration.png",dpi=170)
    summary = {
        "period":args.period,"completed_files":len(tasks),
        "truth_near_tags":totals["count"],"closure_failed":totals["closure_failed"],
        "definition":"Same truth-near selection as kinematic_migration.py; each "
                     "charged track receives its Shapley-average contribution "
                     "to reconstructed minus truth missing-proton polar angle.",
        "notes":"This is exact additive bookkeeping for the chosen attribution "
                "convention, not causal proof or a proton detection efficiency."
    }
    (output/"4_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(f"[OUTPUT] {output} elapsed={time.monotonic()-started:.1f}s",flush=True)


if __name__ == "__main__":
    main()
