#!/usr/bin/env python3
"""Quick truth-near MC diagnostic of electron and missing-proton angle shifts."""
import argparse
from array import array
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

MP2 = .9382720813**2
RAD = 180/np.pi

def part(task):
    path, limit = task
    import ROOT
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    try:
        tree = f.Get("Tags")
        needed = {"truth_tag_match_valid", "truth_missing_M2", "truth_best_index",
                  "truth_best_match_angle", "missing_theta", "truth_missing_theta",
                  "missing_p", "e_theta", "e_phi", "truth_e_theta", "truth_e_phi"}
        missing = needed - {b.GetName() for b in tree.GetListOfBranches()}
        if missing:
            raise RuntimeError(f"Missing branches in {path}: {sorted(missing)}")
        df = ROOT.RDataFrame(tree)
        if limit is not None:
            df = df.Range(0, min(limit, int(tree.GetEntries())))
        ROOT.gInterpreter.Declare("#include <algorithm>")
        df = df.Filter(
            "truth_tag_match_valid>0.5 && "
            f"truth_missing_M2>={MP2-.15:.17g} && truth_missing_M2<{MP2+.15:.17g} && "
            "truth_best_index>=0 && truth_best_match_angle<10 && "
            "std::isfinite(missing_theta) && std::isfinite(truth_missing_theta) && "
            "std::isfinite(missing_p) && missing_p>=0.2 && missing_p<0.5 && "
            "std::isfinite(e_theta) && std::isfinite(e_phi) && "
            "std::isfinite(truth_e_theta) && std::isfinite(truth_e_phi) && "
            "missing_theta>=50 && missing_theta<70")
        df = df.Define("electron_polar_error", "e_theta*57.29577951308232-truth_e_theta")
        df = df.Define("proton_angle_error", "missing_theta-truth_missing_theta")
        df = df.Define("electron_opening_error",
            "std::acos(std::clamp(std::cos(e_theta)*std::cos(truth_e_theta*0.017453292519943295)"
            "+std::sin(e_theta)*std::sin(truth_e_theta*0.017453292519943295)"
            "*std::cos(e_phi-truth_e_phi*0.017453292519943295),-1.0,1.0))"
            "*57.29577951308232")
        out = {}
        handles = {}
        for theta in (50, 60):
            cell = df.Filter(f"missing_theta>={theta} && missing_theta<{theta+10}")
            for key, lo, hi, nb in (("electron_polar_error", -12, 12, 120),
                                    ("electron_opening_error", 0, 20, 100)):
                name = f"h_{theta}_{key}"
                handles[theta,key] = cell.Histo2D(
                    ROOT.RDF.TH2DModel(name, "", nb, lo, hi, 120, -60, 60),
                    key, "proton_angle_error")
        ROOT.RDF.RunGraphs(list(handles.values()))
        for (theta,key), handle in handles.items():
            h = handle.GetValue()
            out[theta,key] = np.array(
                [[h.GetBinContent(i+1,j+1) for j in range(120)]
                 for i in range(h.GetNbinsX())], dtype=np.int64)
        return out
    finally:
        f.Close()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-root", type=Path, required=True)
    ap.add_argument("--period", default="fa18_inb")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--max-tags", type=int, default=10000000,
                    help="0 reads every completed MC part")
    ap.add_argument("--output-root", type=Path, default=Path("output/proton_efficiency"))
    args = ap.parse_args()
    if not 1 <= args.workers <= 10:
        ap.error("workers must be 1–10")
    paths = sorted((args.input_root/args.period/"mc"/"parts").glob("*.root"))
    tasks, remaining = [], args.max_tags
    for p in paths:
        marker = p.with_suffix(".ok.json")
        if not marker.is_file():
            continue
        n = int(json.loads(marker.read_text()).get("selected_tags", 0))
        if n <= 0:
            continue
        limit = min(n,remaining) if args.max_tags else None
        tasks.append((str(p),limit))
        if args.max_tags:
            remaining -= limit
            if remaining == 0:
                break
    if not tasks or (args.max_tags and remaining):
        ap.error("Not enough completed MC tags")
    print(f"[MC] {len(tasks)} files, {min(args.workers,len(tasks))} workers", flush=True)
    total = {}
    with ProcessPoolExecutor(max_workers=min(args.workers,len(tasks))) as pool:
        futures = [pool.submit(part, task) for task in tasks]
        for i, future in enumerate(as_completed(futures),1):
            result = future.result()
            for key, values in result.items():
                total[key] = total.get(key, np.zeros_like(values)) + values
            if i % 10 == 0 or i == len(tasks):
                print(f"[MC] {i}/{len(tasks)} files", flush=True)
    out = args.output_root/args.period/(
        "electron_angle_residuals_full" if not args.max_tags
        else f"electron_angle_residuals_{args.max_tags}")
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2,2,figsize=(12,9),sharey=True,constrained_layout=True)
    for col, theta in enumerate((50,60)):
        for row, (key, lo, hi, label) in enumerate((
            ("electron_polar_error",-12,12,"Electron θrec − θtruth (degrees)"),
            ("electron_opening_error",0,20,"Electron rec–truth opening angle (degrees)"))):
            ax = axes[row,col]
            values = total[theta,key]
            xedges = np.linspace(lo,hi,values.shape[0]+1)
            yedges = np.linspace(-60,60,121)
            image = ax.pcolormesh(xedges,yedges,values.T,
                                  norm=LogNorm(vmin=1),cmap="viridis",shading="flat")
            ax.axhline(0,color="white",lw=.8)
            if row == 0:
                ax.axvline(0,color="white",lw=.8)
            ax.set(xlabel=label,
                   ylabel="Missing-proton θrec − θtruth (degrees)",
                   title=f"Predicted proton θ={theta}–{theta+10}°, p=0.2–0.5 GeV; N={values.sum():,}")
            fig.colorbar(image,ax=ax,label="Tags per cell")
    fig.suptitle("Electron angle error versus predicted proton angle error (truth-near MC)")
    target = out/"1_electron_vs_missing_proton_angle.png"
    fig.savefig(target,dpi=160)
    plt.close(fig)
    print(f"[OUTPUT] {target}",flush=True)

if __name__ == "__main__":
    main()
