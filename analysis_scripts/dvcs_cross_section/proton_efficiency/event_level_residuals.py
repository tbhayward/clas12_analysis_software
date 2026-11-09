#!/usr/bin/env python3
"""MC event-level diagnosis of the proton-tag missing-mass shift.

Uses existing Tags trees. The e pi+ pi- missing mass does not use the detected
proton: its residual is decomposed by measured charged track. Separately, for
angle-matched selected protons, compare reconstructed and generated momentum.
No efficiency correction or energy-loss correction is applied.
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
import numpy as np


PROTON_MASS = 0.9382720813
THETA_EDGES = np.array([0, 10, 20, 30, 40, 50, 60, 70], dtype=float)
P_EDGES = np.array([.2, .5, .7, .9,
                    math.sqrt(1.0 + 1.0 / (4.0 * PROTON_MASS**2))])
M2_DELTA_EDGES = np.linspace(-1.2, 1.2, 241)
P_DELTA_EDGES = np.linspace(-.8, .8, 201)
HISTOGRAMS = {
    "delta_m2": M2_DELTA_EDGES,
    "electron_m2": M2_DELTA_EDGES,
    "pip_m2": M2_DELTA_EDGES,
    "pim_m2": M2_DELTA_EDGES,
    "delta_missing_p": P_DELTA_EDGES,
    "candidate_minus_truth_p": P_DELTA_EDGES,
    "truth_minus_missing_p": P_DELTA_EDGES,
    "candidate_minus_missing_p": P_DELTA_EDGES,
}
PROTON_HISTOGRAMS = ("candidate_minus_truth_p", "truth_minus_missing_p",
                     "candidate_minus_missing_p")
REQUIRED = {
    "missing_M2", "missing_p", "missing_theta", "truth_tag_match_valid",
    "truth_missing_M2", "truth_missing_p", "truth_missing_E", "truth_missing_px",
    "truth_missing_py", "truth_missing_pz", "truth_best_index",
    "truth_best_match_angle", "truth_best_p", "best_selected_rec_index",
    "best_selected_match_angle", "best_selected_p",
}
PARTICLES = ("e", "pip", "pim")
for prefix in PARTICLES:
    REQUIRED.update({f"{prefix}_{part}" for part in ("p", "theta", "phi")})
    REQUIRED.update({f"truth_{prefix}_{part}" for part in ("p", "theta", "phi")})


CPP = r"""
#include <ROOT/RVec.hxx>
#include <cmath>
namespace proton_tag_residual {
struct P4 { double e, x, y, z; };
P4 operator-(P4 a, P4 b) { return {a.e-b.e,a.x-b.x,a.y-b.y,a.z-b.z}; }
double dot(P4 a, P4 b) { return a.e*b.e-a.x*b.x-a.y*b.y-a.z*b.z; }
P4 particle(double p, double theta, double phi, double mass) {
  constexpr double rad=3.14159265358979323846/180.0;
  double s=std::sin(theta*rad);
  return {std::hypot(p,mass),p*s*std::cos(phi*rad),
          p*s*std::sin(phi*rad),p*std::cos(theta*rad)};
}
ROOT::VecOps::RVec<double> split_m2(
    double ep,double et,double ef,double tep,double tet,double tef,
    double pp,double pt,double pf,double tpp,double tpt,double tpf,
    double mp,double mt,double mf,double tmp,double tmt,double tmf,
    double me,double mx,double my,double mz) {
  P4 e=particle(ep,et,ef,0.0005109989461)-particle(tep,tet,tef,0.0005109989461);
  P4 pip=particle(pp,pt,pf,0.13957039)-particle(tpp,tpt,tpf,0.13957039);
  P4 pim=particle(mp,mt,mf,0.13957039)-particle(tmp,tmt,tmf,0.13957039);
  P4 missing{me,mx,my,mz};
  double ce=-2*dot(missing,e)+dot(e,e)+dot(e,pip)+dot(e,pim);
  double cp=-2*dot(missing,pip)+dot(pip,pip)+dot(pip,e)+dot(pip,pim);
  double cm=-2*dot(missing,pim)+dot(pim,pim)+dot(pim,e)+dot(pim,pip);
  return {ce,cp,cm};
}
}
"""
# The C++ signature pairs reconstructed and truth values particle by particle.
SPLIT_ARGS = [value for prefix in PARTICLES
              for value in (f"{prefix}_p", f"{prefix}_theta", f"{prefix}_phi",
                            f"truth_{prefix}_p", f"truth_{prefix}_theta",
                            f"truth_{prefix}_phi")]
SPLIT_ARGS += ["truth_missing_E", "truth_missing_px", "truth_missing_py",
               "truth_missing_pz"]


_CPP_DECLARED = False


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
            raise RuntimeError(f"Missing Tags tree: {filename}")
        branches = {branch.GetName() for branch in tree.GetListOfBranches()}
        missing = REQUIRED - branches
        if missing:
            raise RuntimeError(f"Missing branches in {filename}: {sorted(missing)}")
        if not _CPP_DECLARED:
            ROOT.gInterpreter.Declare(CPP)
            _CPP_DECLARED = True
        frame = ROOT.RDataFrame(tree)
        if stop is not None:
            frame = frame.Range(0, min(stop, int(tree.GetEntries())))
        frame = frame.Filter(
            "truth_tag_match_valid > 0.5 && std::isfinite(truth_missing_M2) && "
            "std::isfinite(missing_M2) && std::isfinite(missing_p) && "
            "std::isfinite(truth_missing_p) && std::isfinite(missing_theta) && "
            "truth_best_index >= 0 && std::isfinite(truth_best_match_angle) && "
            f"truth_best_match_angle < {match_angle:.17g} && "
            f"truth_missing_M2 >= {PROTON_MASS**2-.15:.17g} && "
            f"truth_missing_M2 < {PROTON_MASS**2+.15:.17g} && "
            f"missing_theta >= {THETA_EDGES[0]:.17g} && "
            f"missing_theta < {THETA_EDGES[-1]:.17g} && "
            f"missing_p >= {P_EDGES[0]:.17g} && "
            f"missing_p < {P_EDGES[-1]:.17g}")
        for prefix in PARTICLES:
            for part in ("p", "theta", "phi"):
                frame = frame.Filter(f"std::isfinite({prefix}_{part}) && "
                                     f"std::isfinite(truth_{prefix}_{part})")
        for part in ("E", "px", "py", "pz"):
            frame = frame.Filter(f"std::isfinite(truth_missing_{part})")
        frame = frame.Define("components", "proton_tag_residual::split_m2(" +
                             ",".join(SPLIT_ARGS) + ")")
        frame = (frame.Define("delta_m2", "missing_M2-truth_missing_M2")
                      .Define("electron_m2", "components[0]")
                      .Define("pip_m2", "components[1]")
                      .Define("pim_m2", "components[2]")
                      .Define("delta_missing_p", "missing_p-truth_missing_p")
                      .Define("closure", "delta_m2-electron_m2-pip_m2-pim_m2"))
        # Small differences are expected from stored floating-point momenta.
        closure_bad = frame.Filter("!std::isfinite(closure) || std::abs(closure)>0.005").Count()
        selected_count = frame.Count()
        proton = frame.Filter("best_selected_rec_index >= 0 && "
                              "std::isfinite(best_selected_match_angle) && "
                              f"best_selected_match_angle < {match_angle:.17g} && "
                              "std::isfinite(best_selected_p) && "
                              "std::isfinite(truth_best_p)")
        proton = (proton.Define("candidate_minus_truth_p", "best_selected_p-truth_best_p")
                        .Define("truth_minus_missing_p", "truth_best_p-missing_p")
                        .Define("candidate_minus_missing_p", "best_selected_p-missing_p"))
        proton_count = proton.Count()
        th_bins, p_bins = array("d", THETA_EDGES), array("d", P_EDGES)
        handles = {}
        for name, edges in HISTOGRAMS.items():
            model = ROOT.RDF.TH3DModel(
                "h_" + name, "", len(edges)-1, array("d", edges),
                len(THETA_EDGES)-1, th_bins, len(P_EDGES)-1, p_bins)
            node = proton if name in PROTON_HISTOGRAMS else frame
            handles[name] = node.Histo3D(model, name, "missing_theta", "missing_p")
        ROOT.RDF.RunGraphs([selected_count, proton_count, closure_bad, *handles.values()])

        def values(handle, edges):
            hist = handle.GetValue()
            return np.array([[[int(hist.GetBinContent(i+1,j+1,k+1))
                               for k in range(len(P_EDGES)-1)]
                              for j in range(len(THETA_EDGES)-1)]
                             for i in range(len(edges)-1)], dtype=np.int64)

        return {"selected": int(selected_count.GetValue()),
                "proton": int(proton_count.GetValue()),
                "closure_bad": int(closure_bad.GetValue()),
                "hist": {name: values(handle, HISTOGRAMS[name])
                         for name, handle in handles.items()}}
    finally:
        source.Close()


def completed_tasks(input_root, period, max_tags, match_angle):
    paths = sorted((input_root / period / "mc" / "parts").glob("*.root"))
    tasks, remaining = [], max_tags
    for path in paths:
        marker = path.with_suffix(".ok.json")
        if not marker.is_file():
            continue
        stop = None
        if max_tags:
            count = int(json.loads(marker.read_text()).get("selected_tags", 0))
            if count <= 0:
                continue
            stop = min(count, remaining)
            remaining -= stop
        tasks.append((str(path), stop, match_angle))
        if max_tags and remaining == 0:
            break
    if not tasks or (max_tags and remaining):
        raise RuntimeError("No complete MC input, or fewer tags than requested")
    return tasks


def scan(tasks, workers):
    shape = lambda edges: (len(edges)-1, len(THETA_EDGES)-1, len(P_EDGES)-1)
    total = {"selected": 0, "proton": 0, "closure_bad": 0,
             "hist": {name: np.zeros(shape(edges), dtype=np.int64)
                      for name, edges in HISTOGRAMS.items()}}
    print(f"[MC] files={len(tasks)} workers={min(workers,len(tasks))}", flush=True)
    with ProcessPoolExecutor(max_workers=min(workers, len(tasks))) as pool:
        futures = {pool.submit(read_part, task): task[0] for task in tasks}
        for n, future in enumerate(as_completed(futures), 1):
            try:
                result = future.result()
            except Exception as error:
                raise RuntimeError(f"Failed {futures[future]}: {error}") from error
            for key in ("selected", "proton", "closure_bad"):
                total[key] += result[key]
            for name in HISTOGRAMS:
                total["hist"][name] += result["hist"][name]
            if n % 10 == 0 or n == len(tasks):
                print(f"[MC] {n}/{len(tasks)} files ({100*n/len(tasks):.1f}%) "
                      f"truth-near={total['selected']:,} proton={total['proton']:,}",
                      flush=True)
    return total


def metrics(hist, edges, j, k):
    counts = hist[:, j, k]
    n = int(counts.sum())
    if n == 0:
        return (0, float("nan"), float("nan"))
    centers = (edges[:-1]+edges[1:])/2
    mean = float(np.dot(counts, centers)/n)
    median = float(centers[np.searchsorted(np.cumsum(counts), (n+1)//2)])
    return (n, mean, median)


def write_metrics(total, output):
    with (output / "8_residual_cells.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("theta_low_deg", "theta_high_deg", "p_low_GeV", "p_high_GeV",
                         "observable", "entries_in_histogram_range", "mean", "median"))
        for j in range(len(THETA_EDGES)-1):
            for k in range(len(P_EDGES)-1):
                for name, edges in HISTOGRAMS.items():
                    writer.writerow((THETA_EDGES[j], THETA_EDGES[j+1], P_EDGES[k],
                                     P_EDGES[k+1], name,
                                     *metrics(total["hist"][name], edges, j, k)))


def plot_summary(total, output):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
    names = ("delta_m2", "electron_m2", "pip_m2", "pim_m2")
    titles = ("Total reconstructed − truth M²", "Electron contribution",
              "π+ contribution", "π− contribution")
    for ax, name, title in zip(axes.flat, names, titles):
        values = np.full((len(THETA_EDGES)-1, len(P_EDGES)-1), np.nan)
        for j in range(values.shape[0]):
            for k in range(values.shape[1]):
                n, mean, _ = metrics(total["hist"][name], HISTOGRAMS[name], j, k)
                if n >= 100:
                    values[j,k] = mean
                    ax.text((P_EDGES[k]+P_EDGES[k+1])/2,
                            (THETA_EDGES[j]+THETA_EDGES[j+1])/2, f"{mean:+.2f}",
                            ha="center", va="center", fontsize=8,
                            color="white" if abs(mean)>.45 else "black")
        image = ax.pcolormesh(P_EDGES, THETA_EDGES, values, cmap="coolwarm",
                              vmin=-.8, vmax=.8, shading="flat")
        ax.set(title=title, xlabel="Predicted proton p [GeV]",
               ylabel="Predicted proton θ [degrees]")
        fig.colorbar(image, ax=ax, label="Mean ΔM² [GeV²]")
    fig.suptitle("Truth-near, generated-proton-aligned MC; means of in-range residuals")
    fig.tight_layout()
    fig.savefig(output / "0_track_contribution_summary.png", dpi=150)
    plt.close(fig)


def plot_slices(total, output, theta_index, proton=False):
    names = (PROTON_HISTOGRAMS if proton else
             ("delta_m2", "electron_m2", "pip_m2", "pim_m2"))
    labels = ({"candidate_minus_truth_p": "Reconstructed proton − generated proton",
               "truth_minus_missing_p": "Generated proton − predicted proton",
               "candidate_minus_missing_p": "Reconstructed proton − predicted proton"}
              if proton else
              {"delta_m2": "Total", "electron_m2": "Electron",
               "pip_m2": "π+", "pim_m2": "π−"})
    colors = ("black", "tab:blue", "tab:orange", "tab:green")
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    for k, ax in enumerate(axes.flat):
        for color, name in zip(colors, names):
            hist = total["hist"][name][:, theta_index, k]
            n = int(hist.sum())
            if n:
                ax.stairs(hist/n, HISTOGRAMS[name], color=color,
                          label=f"{labels[name]} (N={n:,})")
        ax.axvline(0, color="gray", linestyle=":")
        ax.set_title(f"Predicted p = {P_EDGES[k]:.2f}–{P_EDGES[k+1]:.2f} GeV")
        ax.set_ylabel("Fraction per bin")
        ax.set_xlabel("Momentum residual [GeV]" if proton else "M² contribution [GeV²]")
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    lo, hi = THETA_EDGES[theta_index:theta_index+2]
    fig.suptitle(("Matched selected proton momentum" if proton else
                  "Missing-mass shift by charged track") + f"; predicted θ={lo:g}–{hi:g}°")
    fig.tight_layout()
    number = theta_index - 3 + (3 if proton else 0)
    stem = "proton_momentum" if proton else "track_contributions"
    fig.savefig(output / f"{number}_{stem}_theta_{lo:g}_{hi:g}.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-root", type=Path, default=Path("output/proton_efficiency"))
    parser.add_argument("--period", default="fa18_inb")
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--max-tags", type=int, default=0,
                        help="0 means all completed MC Tags rows; use a limit for a pilot")
    parser.add_argument("--match-angle", type=float, default=10.0)
    args = parser.parse_args()
    if not 1 <= args.workers <= 10 or args.max_tags < 0 or args.match_angle <= 0:
        parser.error("workers must be 1..10, max-tags nonnegative, match-angle positive")
    output = args.output_root / args.period / (f"event_residuals_{args.max_tags}" if args.max_tags
                                                else "event_residuals_full")
    output.mkdir(parents=True, exist_ok=True)
    tasks = completed_tasks(args.input_root, args.period, args.max_tags, args.match_angle)
    started = time.perf_counter()
    total = scan(tasks, args.workers)
    if total["closure_bad"]:
        print(f"[WARNING] decomposition closure failed for {total['closure_bad']:,} "
              "truth-near tags; check 8_residual_cells.csv before interpretation", flush=True)
    np.savez_compressed(output / "7_residual_histograms.npz",
                        theta_edges=THETA_EDGES, p_edges=P_EDGES,
                        **{f"{name}_edges": edges for name, edges in HISTOGRAMS.items()},
                        **total["hist"])
    plot_summary(total, output)
    for theta_index in (4,5,6):
        plot_slices(total, output, theta_index)
        plot_slices(total, output, theta_index, proton=True)
    write_metrics(total, output)
    (output / "9_summary.json").write_text(json.dumps({
        "period": args.period, "files": len(tasks), "max_tags": args.max_tags,
        "selected_truth_near_aligned": total["selected"],
        "selected_proton_angle_matched": total["proton"],
        "decomposition_closure_bad_gt_0p005_GeV2": total["closure_bad"],
        "truth_near_definition": "truth_missing_M2 within proton mass squared ±0.15 GeV²",
        "proton_match_angle_deg": args.match_angle,
        "note": "Proton momentum comparison uses the selected reconstructed proton "
                "nearest the predicted missing direction; it is not MC::RecMatch proton identity."
    }, indent=2) + "\n")
    print(f"[OUTPUT] {output} elapsed={time.perf_counter()-started:.1f}s", flush=True)


if __name__ == "__main__":
    main()
