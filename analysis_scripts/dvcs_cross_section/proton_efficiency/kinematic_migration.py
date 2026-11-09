#!/usr/bin/env python3
"""MC diagnostic for e pi+ pi- missing-proton migration and tag kinematics.

Uses completed campaign Tags trees. Truth-near means truth missing mass squared
within proton mass squared +/- 0.15 GeV squared, with a generated proton within
10 degrees of the truth missing direction. No correction is calculated.
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


PROTON_MASS = 0.9382720813
THETA = np.array([0, 10, 20, 30, 40, 50, 60, 70], dtype=float)
P = np.array([.2, .5, .7, .9, math.sqrt(1 + 1/(4*PROTON_MASS**2))])
EDGES = {
    "truth_theta": np.linspace(0, 180, 181),
    "angle_difference": np.linspace(-70, 70, 141),
    "opening_angle": np.linspace(0, 90, 91),
    "rec_e_p": np.linspace(2, 10.6, 87),
    "truth_e_p": np.linspace(0, 10.6, 107),
    "electron_delta_p": np.linspace(-3, 3, 151),
    "rec_nu": np.linspace(0, 8.6, 87),
    "rec_q2": np.linspace(0, 10, 101),
    "rec_e_theta_deg": np.linspace(0, 45, 91),
    "delta_m2": np.linspace(-1.2, 1.2, 121),
}
STAGES = ("all", "eb_exists", "selected_exists", "eb_angle_match",
          "selected_angle_match")
REQUIRED = {
    "missing_theta", "missing_p", "missing_M2", "missing_px", "missing_py",
    "missing_pz", "truth_tag_match_valid", "truth_missing_theta",
    "truth_missing_p", "truth_missing_M2", "truth_missing_px",
    "truth_missing_py", "truth_missing_pz", "truth_best_index",
    "truth_best_match_angle", "e_p", "truth_e_p", "e_theta", "beam_energy",
    "Q2", "n_proton_eb", "n_proton_selected", "best_eb_rec_index", "best_eb_match_angle",
    "best_selected_rec_index", "best_selected_match_angle",
}


def read_part(task):
    filename, stop, match_angle = task
    import ROOT

    ROOT.TH1.AddDirectory(False)
    source = ROOT.TFile.Open(filename, "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")
    try:
        tree = source.Get("Tags")
        if not tree:
            raise RuntimeError(f"Missing Tags in {filename}")
        branches = {b.GetName() for b in tree.GetListOfBranches()}
        if REQUIRED - branches:
            raise RuntimeError(f"Missing branches: {sorted(REQUIRED - branches)}")
        frame = ROOT.RDataFrame(tree)
        if stop is not None:
            frame = frame.Range(0, min(stop, int(tree.GetEntries())))
        frame = frame.Filter(
            "truth_tag_match_valid > 0.5 && std::isfinite(truth_missing_M2) && "
            "std::isfinite(truth_missing_theta) && truth_missing_theta >= 0 && "
            "truth_missing_theta < 180 && std::isfinite(truth_missing_p) && "
            "std::isfinite(missing_M2) && std::isfinite(missing_theta) && "
            "std::isfinite(missing_p) && std::isfinite(e_p) && "
            "std::isfinite(truth_e_p) && std::isfinite(e_theta) && "
            "std::isfinite(beam_energy) && std::isfinite(Q2) && "
            "std::isfinite(missing_px) && std::isfinite(missing_py) && "
            "std::isfinite(missing_pz) && std::isfinite(truth_missing_px) && "
            "std::isfinite(truth_missing_py) && std::isfinite(truth_missing_pz) && "
            "truth_best_index >= 0 && std::isfinite(truth_best_match_angle) && "
            f"truth_best_match_angle < {match_angle:.17g} && "
            f"truth_missing_M2 >= {PROTON_MASS**2-.15:.17g} && "
            f"truth_missing_M2 < {PROTON_MASS**2+.15:.17g} && "
            f"missing_theta >= {THETA[0]:.17g} && missing_theta < {THETA[-1]:.17g} && "
            f"missing_p >= {P[0]:.17g} && missing_p < {P[-1]:.17g} && "
            "missing_p > 0 && truth_missing_p > 0")
        frame = (frame.Define("angle_difference", "truth_missing_theta-missing_theta")
                      .Define("electron_delta_p", "e_p-truth_e_p")
                      .Define("rec_nu", "beam_energy-e_p")
                      .Define("rec_e_theta_deg", "e_theta*57.29577951308232")
                      .Define("delta_m2", "missing_M2-truth_missing_M2")
                      .Define("truth_theta", "truth_missing_theta")
                      .Define("rec_e_p", "e_p")
                      .Define("rec_q2", "Q2"))
        cosine = ("(missing_px*truth_missing_px + missing_py*truth_missing_py + "
                  "missing_pz*truth_missing_pz)/(missing_p*truth_missing_p)")
        frame = frame.Define(
            "opening_angle", "std::acos(std::max(-1.0,std::min(1.0," + cosine +
            ")))*57.29577951308232")
        count = frame.Count()
        stages = {"all": frame,
                  "eb_exists": frame.Filter("n_proton_eb > 0"),
                  "selected_exists": frame.Filter("n_proton_selected > 0"),
                  "eb_angle_match": frame.Filter("best_eb_rec_index >= 0 && "
                                                  "std::isfinite(best_eb_match_angle) && "
                                                  f"best_eb_match_angle < {match_angle:.17g}"),
                  "selected_angle_match": frame.Filter(
                      "best_selected_rec_index >= 0 && "
                      "std::isfinite(best_selected_match_angle) && "
                      f"best_selected_match_angle < {match_angle:.17g}")}
        th_bins, p_bins = array("d", THETA), array("d", P)
        counts = {name: node.Histo2D(
            ROOT.RDF.TH2DModel("h_count_" + name, "", len(THETA)-1, th_bins,
                               len(P)-1, p_bins), "missing_theta", "missing_p")
                  for name, node in stages.items()}
        hists = {name: frame.Histo3D(
            ROOT.RDF.TH3DModel("h_" + name, "", len(edges)-1, array("d", edges),
                               len(THETA)-1, th_bins, len(P)-1, p_bins),
            name, "missing_theta", "missing_p") for name, edges in EDGES.items()}
        # Correlations are pooled over momentum within each high-angle band.
        corr = {}
        migration = {}
        for j in (4, 5, 6):
            lo, hi = THETA[j:j+2]
            node = frame.Filter(f"missing_theta >= {lo:.17g} && missing_theta < {hi:.17g}")
            corr[j] = node.Histo2D(
                ROOT.RDF.TH2DModel(f"h_corr_{j}", "", 100, -2., 2., 100, -1.2, 1.2),
                "electron_delta_p", "delta_m2")
        for k in range(len(P)-1):
            lo, hi = P[k:k+2]
            node = frame.Filter(f"missing_p >= {lo:.17g} && missing_p < {hi:.17g}")
            migration[k] = node.Histo2D(
                ROOT.RDF.TH2DModel(f"h_migration_{k}", "", 70, 0., 70., 180, 0., 180.),
                "missing_theta", "truth_missing_theta")
        ROOT.RDF.RunGraphs([count, *counts.values(), *hists.values(),
                            *corr.values(), *migration.values()])

        def h2(handle, nx, ny):
            h = handle.GetValue()
            return np.array([[int(h.GetBinContent(i+1,j+1)) for j in range(ny)]
                             for i in range(nx)], dtype=np.int64)

        def h3(handle, nx):
            h = handle.GetValue()
            return np.array([[[int(h.GetBinContent(i+1,j+1,k+1))
                               for k in range(len(P)-1)] for j in range(len(THETA)-1)]
                             for i in range(nx)], dtype=np.int64)

        return {"selected": int(count.GetValue()),
                "counts": {key: h2(h, len(THETA)-1, len(P)-1)
                           for key, h in counts.items()},
                "hist": {key: h3(h, len(EDGES[key])-1) for key, h in hists.items()},
                "corr": {key: h2(h, 100, 100) for key, h in corr.items()},
                "migration": {key: h2(h, 70, 180) for key, h in migration.items()}}
    finally:
        source.Close()


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


def scan(tasks, workers):
    total = {"selected": 0,
             "counts": {key: np.zeros((len(THETA)-1,len(P)-1), np.int64) for key in STAGES},
             "hist": {key: np.zeros((len(edges)-1,len(THETA)-1,len(P)-1), np.int64)
                      for key, edges in EDGES.items()},
             "corr": {j: np.zeros((100,100), np.int64) for j in (4,5,6)},
             "migration": {k: np.zeros((70,180), np.int64) for k in range(len(P)-1)}}
    print(f"[MC] files={len(tasks)} workers={min(workers,len(tasks))}", flush=True)
    with ProcessPoolExecutor(max_workers=min(workers,len(tasks))) as pool:
        futures = {pool.submit(read_part, task): task[0] for task in tasks}
        for n, future in enumerate(as_completed(futures), 1):
            try:
                part = future.result()
            except Exception as error:
                raise RuntimeError(f"Failed {futures[future]}: {error}") from error
            total["selected"] += part["selected"]
            for group in ("counts","hist","corr","migration"):
                for key in total[group]:
                    total[group][key] += part[group][key]
            if n%10 == 0 or n==len(tasks):
                print(f"[MC] {n}/{len(tasks)} ({100*n/len(tasks):.1f}%) "
                      f"truth-near={total['selected']:,}", flush=True)
    return total


def mean_in_cell(hist, edges, j, k):
    values = hist[:,j,k]
    n = int(values.sum())
    centers = (edges[1:]+edges[:-1])/2
    return (n, float(np.dot(values, centers)/n) if n else float("nan"))


def write_cells(total, output):
    with (output/"6_migration_cells.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("theta_low_deg","theta_high_deg","p_low_GeV","p_high_GeV",
                         "all_tags","eb_exists","selected_exists","eb_angle_match",
                         "selected_angle_match","mean_rec_e_p_GeV","mean_truth_e_p_GeV",
                         "mean_e_p_residual_GeV","mean_virtual_photon_energy_GeV",
                         "mean_rec_Q2_GeV2","mean_rec_e_theta_deg",
                         "mean_truth_minus_rec_proton_theta_deg",
                         "mean_truth_rec_opening_angle_deg","mean_mass_squared_residual_GeV2"))
        for j in range(len(THETA)-1):
            for k in range(len(P)-1):
                values = [mean_in_cell(total["hist"][name], EDGES[name], j, k)[1]
                          for name in ("rec_e_p","truth_e_p","electron_delta_p","rec_nu",
                                       "rec_q2","rec_e_theta_deg","angle_difference",
                                       "opening_angle","delta_m2")]
                writer.writerow((THETA[j],THETA[j+1],P[k],P[k+1],
                                 *(int(total["counts"][key][j,k]) for key in STAGES),
                                 *values))


def heatmap(ax, values, title, cmap="viridis", vmin=None, vmax=None, fmt=".2f"):
    image = ax.pcolormesh(P,THETA,values,cmap=cmap,vmin=vmin,vmax=vmax,shading="flat")
    for j in range(values.shape[0]):
        for k in range(values.shape[1]):
            if np.isfinite(values[j,k]):
                ax.text((P[k]+P[k+1])/2,(THETA[j]+THETA[j+1])/2,
                        format(values[j,k],fmt),ha="center",va="center",fontsize=8)
    ax.set(title=title,xlabel="Predicted proton p [GeV]",
           ylabel="Predicted proton θ [degrees]")
    return image


def plot_kinematics(total, output):
    panels = (("rec_e_p","Reconstructed electron p [GeV]"),
              ("rec_nu","Virtual-photon energy ν = Ebeam − p_e [GeV]"),
              ("electron_delta_p","Reconstructed − truth electron p [GeV]"),
              ("rec_q2","Reconstructed Q² [GeV²]"))
    fig, axes = plt.subplots(2,2,figsize=(13,9),sharex=True,sharey=True)
    for ax,(key,title) in zip(axes.flat,panels):
        values = np.full((len(THETA)-1,len(P)-1),np.nan)
        for j in range(values.shape[0]):
            for k in range(values.shape[1]):
                n, mean = mean_in_cell(total["hist"][key],EDGES[key],j,k)
                if n >= 100:
                    values[j,k] = mean
        image = heatmap(ax,values,title)
        fig.colorbar(image,ax=ax,label="Mean in plotted range")
    fig.suptitle("Electron kinematics in reconstructed missing-proton cells")
    fig.tight_layout()
    fig.savefig(output/"1_electron_kinematics.png",dpi=150)
    plt.close(fig)


def plot_migration(total, output):
    fig,axes=plt.subplots(2,2,figsize=(12,10),sharex=True,sharey=True)
    for k,ax in enumerate(axes.flat):
        counts=total["migration"][k]
        image=ax.pcolormesh(np.linspace(0,70,71),np.linspace(0,180,181),counts.T,
                            norm=LogNorm(vmin=1),cmap="viridis",shading="flat")
        ax.plot([0,70],[0,70],color="white",linestyle="--",linewidth=1)
        ax.set_ylim(0,90)
        ax.set(title=f"Predicted p={P[k]:.2f}–{P[k+1]:.2f} GeV; N={counts.sum():,}",
               xlabel="Reconstructed predicted proton θ [degrees]",
               ylabel="Truth missing-proton θ [degrees]")
        fig.colorbar(image,ax=ax,label="Tags / bin")
    fig.suptitle("Truth versus reconstructed predicted proton angle")
    fig.tight_layout()
    fig.savefig(output/"0_proton_angle_migration.png",dpi=150)
    plt.close(fig)


def plot_electron_spectra(total, output):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, j in zip(axes, (4, 5, 6)):
        for key, label, color in (("rec_e_p", "Reconstructed electron", "tab:blue"),
                                  ("truth_e_p", "Matched truth electron", "tab:orange")):
            counts = total["hist"][key][:, j, :].sum(axis=1)
            n = int(counts.sum())
            if n:
                ax.stairs(counts/n, EDGES[key], color=color,
                          label=f"{label} (N={n:,})")
        ax.set(title=f"Predicted proton θ={THETA[j]:g}–{THETA[j+1]:g}°",
               xlabel="Electron momentum [GeV]", ylabel="Fraction / momentum bin")
        ax.set_xlim(0, 10.6)
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle("Electron momentum; lower p_e corresponds to higher virtual-photon energy")
    fig.tight_layout()
    fig.savefig(output/"2_electron_momentum_spectra.png", dpi=150)
    plt.close(fig)


def plot_electron_mass(total,output):
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),sharex=True,sharey=True)
    for ax,j in zip(axes,(4,5,6)):
        counts=total["corr"][j]
        image=ax.pcolormesh(np.linspace(-2,2,101),np.linspace(-1.2,1.2,101),
                            counts.T,norm=LogNorm(vmin=1),cmap="viridis",shading="flat")
        ax.axhline(0,color="white",linewidth=.7)
        ax.axvline(0,color="white",linewidth=.7)
        ax.set(title=f"Predicted θ={THETA[j]:g}–{THETA[j+1]:g}°; N={counts.sum():,}",
               xlabel="Electron p_rec − p_truth [GeV]",
               ylabel="Missing M²_rec − M²_truth [GeV²]")
        fig.colorbar(image,ax=ax,label="Tags / bin")
    fig.suptitle("Electron momentum residual versus missing-mass shift")
    fig.tight_layout()
    fig.savefig(output/"3_electron_residual_vs_mass.png",dpi=150)
    plt.close(fig)


def plot_cutflow(total,output):
    denominator=total["counts"]["all"]
    panels=(("eb_exists","Event Builder proton present"),
            ("selected_exists","Selected proton present"),
            ("eb_angle_match","EB proton within 10° of prediction"),
            ("selected_angle_match","Selected proton within 10° of prediction"))
    fig,axes=plt.subplots(2,2,figsize=(13,9),sharex=True,sharey=True)
    for ax,(key,title) in zip(axes.flat,panels):
        values=np.divide(total["counts"][key],denominator,
                         out=np.full(denominator.shape,np.nan),where=denominator>=100)
        image=heatmap(ax,values,title,vmin=0,vmax=1)
        fig.colorbar(image,ax=ax,label="Fraction of truth-near aligned tags")
    fig.suptitle("Reconstructed proton presence and matching; these are not efficiencies")
    fig.tight_layout()
    fig.savefig(output/"4_proton_candidate_cutflow.png",dpi=150)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root",required=True,type=Path)
    parser.add_argument("--output-root",type=Path,default=Path("output/proton_efficiency"))
    parser.add_argument("--period",default="fa18_inb")
    parser.add_argument("--workers",type=int,default=10)
    parser.add_argument("--max-tags",type=int,default=0,
                        help="0 scans every completed MC part")
    parser.add_argument("--match-angle",type=float,default=10.)
    args=parser.parse_args()
    if not 1<=args.workers<=10 or args.max_tags<0 or args.match_angle<=0:
        parser.error("workers must be 1..10; max-tags nonnegative; match-angle positive")
    tasks=completed_tasks(args.input_root,args.period,args.max_tags,args.match_angle)
    output=args.output_root/args.period/(f"migration_{args.max_tags}" if args.max_tags
                                         else "migration_full")
    output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    total=scan(tasks,args.workers)
    if not np.array_equal(total["counts"]["all"],
                          total["hist"]["truth_theta"].sum(axis=0)):
        raise RuntimeError("Truth-theta histogram does not match tag counts")
    np.savez_compressed(output/"5_migration_histograms.npz",theta_edges=THETA,p_edges=P,
                        **{f"{key}_edges":value for key,value in EDGES.items()},
                        **total["hist"],
                        **{f"count_{key}":value for key,value in total["counts"].items()},
                        **{f"corr_{key}":value for key,value in total["corr"].items()},
                        **{f"migration_{key}":value for key,value in total["migration"].items()})
    plot_migration(total,output)
    plot_kinematics(total,output)
    plot_electron_spectra(total,output)
    plot_electron_mass(total,output)
    plot_cutflow(total,output)
    write_cells(total,output)
    (output/"7_summary.json").write_text(json.dumps({
        "period":args.period,"completed_files":len(tasks),"requested_tags":args.max_tags,
        "truth_near_aligned_tags":total["selected"],"match_angle_deg":args.match_angle,
        "counts_by_stage":{key:int(value.sum()) for key,value in total["counts"].items()},
        "definition":"Truth charged tracks matched via MC::RecMatch, truth missing "
                     "M2 within proton mass squared +/-0.15 GeV2, and a generated "
                     "proton within match_angle_deg of truth missing direction",
        "notes":"Cells use reconstructed predicted proton theta and momentum. "
                "Ebeam-e_p is virtual-photon energy transfer, not a measured final-state photon. "
                "Histogram means exclude underflow/overflow; candidate rates are not efficiencies."
    },indent=2)+"\n")
    print(f"[OUTPUT] {output} elapsed={time.perf_counter()-started:.1f}s",flush=True)


if __name__=="__main__":
    main()
