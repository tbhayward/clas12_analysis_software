#!/usr/bin/env python3
"""Read-only |t|-dependent exclusivity diagnostics; run from external_scripts/.
Requires numpy, matplotlib, uproot. No caches, no PDFs.
Data cut-survival is RAW (signal+background), not a DVCS efficiency.
"""
import argparse
import csv
import json
import os
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import plot_exclusivity_data_dvcs_pi0_mc as original

VARS = ("Delta_phi", "theta", "theta_gamma_gamma", "pTmiss",
        "Emiss2", "Mx2_1", "Mx2")
ALIASES = {"theta": ("theta", "theta2", "theta_2"),
           "Mx2_1": ("Mx2_1",), "Mx2": ("Mx2", "Mx2_epg", "Mx2_eg", "Mx2_epgamma")}
COLORS = {"data": "black", "mc": "tab:blue", "pi0": "tab:orange"}

def get_edges(path):
    # Exact |t| edges used by the production 4D binning CSV, collapsed
    # across xB, Q2 and phi for this first resolution diagnostic.
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        vals = sorted({float(row[k]) for row in r
                       for k in ("t_abs_min", "t_abs_max") if row.get(k)})
    if len(vals) < 2:
        raise RuntimeError("No production |t| bin edges found")
    return np.array(vals)

def resolve(path):
    with uproot.open(f"{path}:PhysicsEvents") as tree:
        keys = set(tree.keys())
    needed = ["detector1", "detector2", "t1", "open_angle_ep2",
              "pTmiss", "e_phi", "p1_phi", "p2_phi", "p1_p"]
    out = {}
    for v in VARS:
        options = ALIASES.get(v, (v,))
        match = next((k for k in options if k in keys), None)
        if match is not None:
            out[v] = match
    for name in needed:
        if name not in keys:
            raise RuntimeError(f"{path}: missing mandatory branch {name}")
    # t2 may or may not represent the independent photon-derived t.
    for candidate in ("t2",):
        if candidate in keys:
            out["t_photon_candidate"] = candidate
    return sorted(set(needed + list(out.values()))), out

def bounds(stat):
    if "cut_low" in stat and "cut_high" in stat:
        return float(stat["cut_low"]), float(stat["cut_high"])
    mu, sigma = float(stat["mean"]), float(stat["std"])
    return mu-3*sigma, mu+3*sigma

def worker(job):
    period, topology, kind, path, cutmap, edges, step = job
    branches, names = resolve(path)
    nt = len(edges)-1
    store = {v: [[] for _ in range(nt)] for v in list(VARS)+["delta_t", "p1_p"]}
    before = np.zeros(nt, np.int64)
    after = np.zeros(nt, np.int64)
    single = {v: np.zeros(nt, np.int64) for v in cutmap}
    # Store only preselected events, in float32, and avoid retaining full ROOT arrays.
    for a in uproot.iterate(f"{path}:PhysicsEvents", expressions=branches,
                            step_size=step, library="np"):
        det1, det2 = a["detector1"], a["detector2"]
        t = -np.asarray(a["t1"])
        mask = ((det1 == topology.detector1) & (det2 == topology.detector2)
                & np.isfinite(t) & (t < 1.0)
                & (np.asarray(a["open_angle_ep2"]) > 5.0))
        # Match the explicit FD-sector separation used by the production
        # exclusivity-template script, not a generic phi separation.
        es = original.fd_sector_from_phi(np.asarray(a["e_phi"], dtype=float))
        ps = original.fd_sector_from_phi(np.asarray(a["p1_phi"], dtype=float))
        gs = original.fd_sector_from_phi(np.asarray(a["p2_phi"], dtype=float))
        mask &= es >= 1
        if topology.detector1 == 1:
            mask &= (ps >= 1) & (ps != es)
        if topology.detector2 == 1:
            mask &= (gs >= 1) & (gs != es)
        if topology.detector1 == 1 and topology.detector2 == 1:
            mask &= ps != gs
        idx = np.searchsorted(edges, t, side="right")-1
        mask &= (idx >= 0) & (idx < nt)
        if not mask.any():
            continue
        allpass = mask.copy()
        for v, stat in cutmap.items():
            branch = names.get(v)
            if branch is None:
                raise RuntimeError(f"{path}: cut variable {v} not found")
            low, high = bounds(stat)
            values = np.asarray(a[branch])
            good = np.isfinite(values) & (values >= low) & (values <= high)
            single[v] += np.bincount(idx[mask & good], minlength=nt)[:nt]
            allpass &= good
        before += np.bincount(idx[mask], minlength=nt)[:nt]
        after += np.bincount(idx[allpass], minlength=nt)[:nt]
        delta = (np.asarray(a[names["t_photon_candidate"]])-a["t1"]
                 if "t_photon_candidate" in names else None)
        for i in np.unique(idx[mask]):
            sel = mask & (idx == i)
            for v in VARS:
                if v in names:
                    store[v][i].append(np.asarray(a[names[v]][sel], dtype=np.float32))
            store["p1_p"][i].append(np.asarray(a["p1_p"][sel], dtype=np.float32))
            if delta is not None:
                store["delta_t"][i].append(np.asarray(delta[sel], dtype=np.float32))
    arrays = {}
    for v, bins in store.items():
        arrays[v] = [np.concatenate(chunks) if chunks else np.array([], dtype=np.float32)
                     for chunks in bins]
    return period, topology.key, kind, before, after, single, arrays

def summarize(x):
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.nan, np.nan, np.nan, len(x)
    q16, q50, q84 = np.quantile(x, [0.16, 0.50, 0.84])
    return float(q50), float((q84-q16)/2), float(np.std(x)), len(x)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--periods", nargs="+", default=["fa18_inb", "sp18_out"])
    ap.add_argument("--topologies", nargs="+", default=["CD_FD", "CD_FT", "FD_FD"])
    ap.add_argument("--include-pi0", action="store_true")
    ap.add_argument("--step-size", default="100 MB")
    ap.add_argument("--binning", default=str(ROOT/"imports/all_bin_v3.csv"))
    ap.add_argument("--cuts", default=str(ROOT/"output/jsons/combined_cuts.json"))
    ap.add_argument("--output", default=str(ROOT/"output/t_dependent_exclusivity_diagnostics"))
    args = ap.parse_args()
    edges = get_edges(args.binning)
    with open(args.cuts) as f:
        cuts = json.load(f)
    periods = {p.key:p for p in original.PERIODS}
    topologies = {t.key:t for t in original.TOPOLOGIES}
    jobs = []
    for pk in args.periods:
        p = periods[pk]
        for tk in args.topologies:
            t = topologies[tk]
            key = original.cut_json_key("dvcs", p, t)
            if key not in cuts:
                raise KeyError(f"Missing {key} in {args.cuts}")
            for kind, path in [("data",p.data_file), ("mc",p.dvcs_mc_file)] + (
                    [("pi0",p.pi0_as_dvcs_mc_file)] if args.include_pi0 else []):
                if not Path(path).is_file():
                    raise FileNotFoundError(path)
                sample = "mc" if kind != "data" else "data"
                jobs.append((pk,t,kind,path,cuts[key][sample],edges,args.step_size))
    out = Path(args.output)
    out.mkdir(parents=True,exist_ok=True)
    results = {}
    with ProcessPoolExecutor(max_workers=max(1,min(8,args.workers,len(jobs)))) as pool:
        futures = [pool.submit(worker,job) for job in jobs]
        for future in as_completed(futures):
            period, topo, kind, before, after, single, arrays = future.result()
            results[(period,topo,kind)] = (before,after,single,arrays)
            print(f"[done] {period} {topo} {kind}: {before.sum():,} preselected",flush=True)
    rows, widths = [], []
    for (period,topo,kind),(before,after,single,arrays) in sorted(results.items()):
        for i in range(len(edges)-1):
            eff = after[i]/before[i] if before[i] else np.nan
            row = dict(period=period,topology=topo,sample=kind,
                       t_min=edges[i],t_max=edges[i+1],n_before=before[i],
                       n_after=after[i],raw_survival=eff,
                       raw_survival_binomial_error=np.sqrt(eff*(1-eff)/before[i])
                       if before[i] else np.nan)
            for v, counts in single.items():
                row[f"single_{v}"] = counts[i]/before[i] if before[i] else np.nan
            rows.append(row)
            for v, bins in arrays.items():
                median, robust, rms, n = summarize(bins[i])
                widths.append(dict(period=period,topology=topo,sample=kind,
                                   t_min=edges[i],t_max=edges[i+1],variable=v,
                                   median=median,robust_sigma=robust,rms=rms,n=n))
    def write_csv(path, records):
        fields = sorted(set().union(*(r.keys() for r in records)))
        with open(path,"w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)
    write_csv(out/"cut_survival.csv",rows)
    write_csv(out/"resolution_widths.csv",widths)
    for pk in args.periods:
        for tk in args.topologies:
            fig,ax=plt.subplots(figsize=(8,5))
            for kind in ("data","mc","pi0"):
                entry=results.get((pk,tk,kind))
                if entry is None: continue
                pre,post,_,_=entry
                eff=np.divide(post,pre,out=np.full(len(pre),np.nan),where=pre>0)
                ax.errorbar((edges[:-1]+edges[1:])/2,eff,
                            yerr=np.sqrt(np.maximum(0,eff*(1-eff))/np.maximum(pre,1)),
                            fmt="o-",label=kind,color=COLORS[kind])
            ax.set(xlabel=r"$|t|$ (GeV$^2$)",ylabel="Raw exclusivity-cut survival",
                   title=f"{pk} / {tk}")
            ax.set_ylim(0,1.05); ax.grid(alpha=.2);ax.legend()
            fig.tight_layout()
            fig.savefig(out/f"survival_{pk}_{tk}.png",dpi=160)
            plt.close(fig)
            for v in list(VARS)+["delta_t","p1_p"]:
                fig,ax=plt.subplots(figsize=(8,5))
                for kind in ("data","mc","pi0"):
                    entry=results.get((pk,tk,kind))
                    if entry is None: continue
                    arr=entry[3][v]
                    widths_=np.array([summarize(x)[1] for x in arr])
                    ax.plot((edges[:-1]+edges[1:])/2,widths_,"o-",
                            label=kind,color=COLORS[kind])
                ax.set(xlabel=r"$|t|$ (GeV$^2$)",ylabel=f"Robust width: {v}",
                       title=f"{pk} / {tk} / {v}")
                ax.grid(alpha=.2);ax.legend();fig.tight_layout()
                fig.savefig(out/f"width_{pk}_{tk}_{v}.png",dpi=160)
                plt.close(fig)
    with open(out/"README.txt","w") as f:
        f.write("Diagnostic only; no production cuts or cross sections modified.\n"
                "t bins: unique exact |t| edges from imports/all_bin_v3.csv, pooled across xB/Q2.\n"
                "raw_survival (DATA) is NOT signal efficiency: pi0/background contamination remains.\n"
                "MC raw_survival is reconstructed signal survival conditional on preselection.\n"
                "Do NOT use data/MC raw survival ratio as a cross-section correction.\n"
                "Widths are central (q84-q16)/2 of preselected distributions, not fitted signal widths.\n"
                "delta_t only if ROOT branch t2 exists; interpretation must be verified.\n"
                "The global preselection matches the standalone template script's minimal selection;\n"
                "additional production fiducial/global cuts must be audited before claiming closure.\n"
                "No generated-truth migration matrix or data template-fit efficiency is included.\n")
    print(f"[output] {out}")
if __name__=="__main__":
    main()
