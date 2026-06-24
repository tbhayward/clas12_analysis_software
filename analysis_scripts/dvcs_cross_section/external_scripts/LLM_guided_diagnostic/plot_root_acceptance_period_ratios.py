#!/usr/bin/env python3

import os
import math
import argparse
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt

try:
    import uproot
except ImportError:
    raise SystemExit(
        "This script requires uproot. Install with: pip install uproot awkward"
    )
#endtry


OUTPUT_DIR_DEFAULT = "output/individual_acceptance_comparison"
MAX_WORKERS_DEFAULT = 5


FILES = {
    "Fa18 Inb": {
        "gen": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_inb_50nA_10604MeV.root",
        "rec": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV.root",
    },
    "Fa18 Out": {
        "gen": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_50nA_10604MeV.root",
        "rec": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV.root",
    },
    "Sp18 Inb": {
        "gen": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_50nA_10594MeV.root",
        "rec": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_50nA_10594MeV.root",
    },
    "Sp18 Out": {
        "gen": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_out_45nA_10594MeV.root",
        "rec": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_45nA_10594MeV.root",
    },
}


PERIOD_RATIOS = [
    ("Sp18 Inb", "Fa18 Inb", "Sp18 Inb / Fa18 Inb acceptance"),
    ("Sp18 Out", "Fa18 Out", "Sp18 Out / Fa18 Out acceptance"),
]


@dataclass
class VariableConfig:
    key: str
    branch: str
    title: str
    xlabel: str
    nbins: int
    xmin: float
    xmax: float
    use_abs: bool = False
    wrap_phi_deg: bool = False
    radians_to_degrees: bool = False


VARIABLES = [
    VariableConfig(
        key="xB",
        branch="x",
        title=r"$x_B$",
        xlabel=r"$x_B$",
        nbins=24,
        xmin=0.05,
        xmax=0.65,
    ),
    VariableConfig(
        key="Q2",
        branch="Q2",
        title=r"$Q^2$",
        xlabel=r"$Q^2$ (GeV$^2$)",
        nbins=25,
        xmin=1.0,
        xmax=6.0,
    ),
    VariableConfig(
        key="t_abs",
        branch="t",
        title=r"$|t|$",
        xlabel=r"$|t|$ (GeV$^2$)",
        nbins=24,
        xmin=0.0,
        xmax=1.2,
        use_abs=True,
    ),
    VariableConfig(
        key="phi",
        branch="phi2",
        title=r"$\phi$",
        xlabel=r"$\phi$ (deg)",
        nbins=24,
        xmin=0.0,
        xmax=360.0,
        wrap_phi_deg=True,
    ),
    VariableConfig(
        key="theta_e",
        branch="e_theta",
        title=r"$\theta_e$",
        xlabel=r"$\theta_e$ (deg)",
        nbins=24,
        xmin=5.0,
        xmax=35.0,
        radians_to_degrees=True,
    ),
    VariableConfig(
        key="theta_p",
        branch="p1_theta",
        title=r"$\theta_p$",
        xlabel=r"$\theta_p$ (deg)",
        nbins=24,
        xmin=15.0,
        xmax=75.0,
        radians_to_degrees=True,
    ),
    VariableConfig(
        key="theta_gamma",
        branch="p2_theta",
        title=r"$\theta_\gamma$",
        xlabel=r"$\theta_\gamma$ (deg)",
        nbins=24,
        xmin=0.0,
        xmax=35.0,
        radians_to_degrees=True,
    ),
    VariableConfig(
        key="W",
        branch="W",
        title=r"$W$",
        xlabel=r"$W$ (GeV)",
        nbins=24,
        xmin=1.8,
        xmax=4.5,
    ),
    VariableConfig(
        key="y",
        branch="y",
        title=r"$y$",
        xlabel=r"$y$",
        nbins=24,
        xmin=0.0,
        xmax=0.9,
    ),
    VariableConfig(
        key="Mx2",
        branch="Mx2",
        title=r"$M_X^2$",
        xlabel=r"$M_X^2$ (GeV$^2$)",
        nbins=24,
        xmin=-1.0,
        xmax=2.0,
    ),
    VariableConfig(
        key="Emiss2",
        branch="Emiss2",
        title=r"$E_{\mathrm{miss}}^2$",
        xlabel=r"$E_{\mathrm{miss}}^2$ (GeV$^2$)",
        nbins=24,
        xmin=-2.0,
        xmax=4.0,
    ),
    VariableConfig(
        key="pTmiss",
        branch="pTmiss",
        title=r"$p_T^{\mathrm{miss}}$",
        xlabel=r"$p_T^{\mathrm{miss}}$ (GeV)",
        nbins=24,
        xmin=0.0,
        xmax=1.0,
    ),
]


def find_tree_name(path: str, requested_tree: Optional[str] = None) -> str:
    with uproot.open(path) as f:
        if requested_tree is not None:
            if requested_tree in f:
                return requested_tree
            #endif

            if requested_tree + ";1" in f:
                return requested_tree + ";1"
            #endif

            raise RuntimeError(f"Requested tree '{requested_tree}' not found in {path}")
        #endif

        for key, obj in f.items():
            class_name = obj.classname

            if "TTree" in class_name:
                return key
            #endif
        #endfor
    #endwith

    raise RuntimeError(f"No TTree found in {path}")


def auto_radians_to_degrees_if_needed(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    finite = x[np.isfinite(x)]

    if finite.size == 0:
        return x
    #endif

    max_abs = np.nanmax(np.abs(finite))

    if max_abs <= 2.0 * math.pi + 0.1:
        return x * 180.0 / math.pi
    #endif

    return x


def transform_values(values: np.ndarray, cfg: VariableConfig) -> np.ndarray:
    x = np.asarray(values, dtype=float)

    if cfg.radians_to_degrees:
        x = auto_radians_to_degrees_if_needed(x)
    #endif

    if cfg.use_abs:
        x = np.abs(x)
    #endif

    if cfg.wrap_phi_deg:
        x = auto_radians_to_degrees_if_needed(x)
        x = np.mod(x + 360.0, 360.0)
    #endif

    return x


def hist_file_task(
    period: str,
    sample_type: str,
    path: str,
    tree_name: Optional[str],
    variable_configs: List[VariableConfig],
    step_size: str,
) -> Dict:
    if not os.path.exists(path):
        raise RuntimeError(f"Missing file for {period} {sample_type}: {path}")
    #endif

    actual_tree = find_tree_name(path, tree_name)

    out = {
        "period": period,
        "sample_type": sample_type,
        "path": path,
        "tree": actual_tree,
        "n_entries_tree": 0,
        "variables": {},
    }

    with uproot.open(path) as f:
        tree = f[actual_tree]
        available_branches = set(tree.keys())
        out["n_entries_tree"] = int(tree.num_entries)

        usable_configs = [
            cfg for cfg in variable_configs
            if cfg.branch in available_branches
        ]

        missing = [
            cfg.branch for cfg in variable_configs
            if cfg.branch not in available_branches
        ]

        if missing:
            out["missing_branches"] = sorted(set(missing))
        else:
            out["missing_branches"] = []
        #endif

        for cfg in usable_configs:
            edges = np.linspace(cfg.xmin, cfg.xmax, cfg.nbins + 1)

            out["variables"][cfg.key] = {
                "branch": cfg.branch,
                "edges": edges,
                "counts": np.zeros(cfg.nbins, dtype=float),
                "n_finite": 0,
                "n_in_range": 0,
                "finite_min": np.nan,
                "finite_max": np.nan,
            }
        #endfor

        if not usable_configs:
            return out
        #endif

        read_branches = sorted(set(cfg.branch for cfg in usable_configs))

        for arrays in tree.iterate(
            expressions=read_branches,
            library="np",
            step_size=step_size,
        ):
            for cfg in usable_configs:
                vals = arrays[cfg.branch]
                vals = transform_values(vals, cfg)

                finite_mask = np.isfinite(vals)
                vals_finite = vals[finite_mask]

                if vals_finite.size == 0:
                    continue
                #endif

                edges = out["variables"][cfg.key]["edges"]
                counts, _ = np.histogram(vals_finite, bins=edges)

                in_range = (
                    (vals_finite >= cfg.xmin)
                    & (vals_finite < cfg.xmax)
                )

                old_min = out["variables"][cfg.key]["finite_min"]
                old_max = out["variables"][cfg.key]["finite_max"]

                batch_min = float(np.nanmin(vals_finite))
                batch_max = float(np.nanmax(vals_finite))

                if not np.isfinite(old_min):
                    new_min = batch_min
                else:
                    new_min = min(old_min, batch_min)
                #endif

                if not np.isfinite(old_max):
                    new_max = batch_max
                else:
                    new_max = max(old_max, batch_max)
                #endif

                out["variables"][cfg.key]["counts"] += counts.astype(float)
                out["variables"][cfg.key]["n_finite"] += int(vals_finite.size)
                out["variables"][cfg.key]["n_in_range"] += int(np.sum(in_range))
                out["variables"][cfg.key]["finite_min"] = new_min
                out["variables"][cfg.key]["finite_max"] = new_max
            #endfor
        #endfor
    #endwith

    return out


def compute_acceptance(
    gen_counts: np.ndarray,
    rec_counts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    gen = np.asarray(gen_counts, dtype=float)
    rec = np.asarray(rec_counts, dtype=float)

    acc = np.full_like(gen, np.nan, dtype=float)
    acc_err = np.full_like(gen, np.nan, dtype=float)

    valid = gen > 0.0
    acc[valid] = rec[valid] / gen[valid]

    good = valid & (acc >= 0.0)
    acc_clipped = np.clip(acc[good], 0.0, 1.0)

    acc_err[good] = np.sqrt(acc_clipped * (1.0 - acc_clipped) / gen[good])

    above_one = valid & (acc > 1.0)
    rel2 = np.zeros_like(gen, dtype=float)

    rec_pos = rec > 0.0

    rel2[above_one & rec_pos] += 1.0 / rec[above_one & rec_pos]
    rel2[above_one] += 1.0 / gen[above_one]

    acc_err[above_one] = acc[above_one] * np.sqrt(rel2[above_one])

    return acc, acc_err


def ratio_with_error(
    num: np.ndarray,
    den: np.ndarray,
    num_err: np.ndarray,
    den_err: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    ratio = np.full_like(num, np.nan, dtype=float)
    ratio_err = np.full_like(num, np.nan, dtype=float)

    valid = (
        np.isfinite(num)
        & np.isfinite(den)
        & (num > 0.0)
        & (den > 0.0)
    )

    ratio[valid] = num[valid] / den[valid]

    rel2 = np.zeros_like(num, dtype=float)

    num_good = valid & np.isfinite(num_err) & (num_err > 0.0)
    den_good = valid & np.isfinite(den_err) & (den_err > 0.0)

    rel2[num_good] += (num_err[num_good] / num[num_good]) ** 2
    rel2[den_good] += (den_err[den_good] / den[den_good]) ** 2

    ratio_err[valid] = ratio[valid] * np.sqrt(rel2[valid])

    return ratio, ratio_err


def integrated_acceptance_from_counts(
    gen_counts: np.ndarray,
    rec_counts: np.ndarray,
) -> Tuple[float, float, float, float]:
    gen_total = float(np.nansum(gen_counts))
    rec_total = float(np.nansum(rec_counts))

    if gen_total <= 0.0:
        return gen_total, rec_total, np.nan, np.nan
    #endif

    acc = rec_total / gen_total

    acc_clipped = min(max(acc, 0.0), 1.0)
    acc_err = math.sqrt(acc_clipped * (1.0 - acc_clipped) / gen_total)

    if acc > 1.0:
        rel2 = 0.0

        if rec_total > 0.0:
            rel2 += 1.0 / rec_total
        #endif

        rel2 += 1.0 / gen_total
        acc_err = acc * math.sqrt(rel2)
    #endif

    return gen_total, rec_total, acc, acc_err


def make_variable_plot(
    cfg: VariableConfig,
    acceptances: Dict,
    output_dir: str,
) -> List[Dict]:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

    rows = []

    for ipad, (num_period, den_period, panel_title) in enumerate(PERIOD_RATIOS):
        ax = axes[ipad]

        gen_num = acceptances[num_period][cfg.key]["gen_counts"]
        rec_num = acceptances[num_period][cfg.key]["rec_counts"]
        acc_num = acceptances[num_period][cfg.key]["acceptance"]
        acc_num_err = acceptances[num_period][cfg.key]["acceptance_err"]

        gen_den = acceptances[den_period][cfg.key]["gen_counts"]
        rec_den = acceptances[den_period][cfg.key]["rec_counts"]
        acc_den = acceptances[den_period][cfg.key]["acceptance"]
        acc_den_err = acceptances[den_period][cfg.key]["acceptance_err"]

        ratio, ratio_err = ratio_with_error(
            acc_num,
            acc_den,
            acc_num_err,
            acc_den_err,
        )

        edges = acceptances[num_period][cfg.key]["edges"]
        centers = 0.5 * (edges[:-1] + edges[1:])

        _, _, integrated_num, integrated_num_err = integrated_acceptance_from_counts(gen_num, rec_num)
        _, _, integrated_den, integrated_den_err = integrated_acceptance_from_counts(gen_den, rec_den)

        integrated_ratio_arr, integrated_ratio_err_arr = ratio_with_error(
            np.array([integrated_num]),
            np.array([integrated_den]),
            np.array([integrated_num_err]),
            np.array([integrated_den_err]),
        )

        integrated_ratio = integrated_ratio_arr[0]
        integrated_ratio_err = integrated_ratio_err_arr[0]

        ax.errorbar(
            centers,
            ratio,
            yerr=ratio_err,
            marker="o",
            linestyle="-",
            linewidth=1.6,
            markersize=5,
            capsize=2.5,
            label=panel_title,
        )

        if np.isfinite(integrated_ratio):
            ax.axhline(
                integrated_ratio,
                linestyle="--",
                linewidth=1.0,
                alpha=0.8,
                label=f"integrated = {integrated_ratio:.3f}",
            )
        #endif

        ax.axhline(1.0, linestyle=":", linewidth=1.2)
        ax.set_title(panel_title, fontsize=13)
        ax.set_xlabel(cfg.xlabel, fontsize=12)
        ax.set_ylabel("acceptance ratio", fontsize=12)
        ax.set_ylim(0.0, 2.5)
        ax.grid(True, alpha=0.35)
        ax.legend(loc="best", fontsize=9)

        for ibin in range(len(centers)):
            rows.append({
                "variable": cfg.key,
                "branch": cfg.branch,
                "bin_index": ibin,
                "bin_low": edges[ibin],
                "bin_high": edges[ibin + 1],
                "bin_center": centers[ibin],
                "numerator_period": num_period,
                "denominator_period": den_period,
                "ratio_label": panel_title,
                "gen_num": gen_num[ibin],
                "rec_num": rec_num[ibin],
                "acceptance_num": acc_num[ibin],
                "acceptance_num_err": acc_num_err[ibin],
                "gen_den": gen_den[ibin],
                "rec_den": rec_den[ibin],
                "acceptance_den": acc_den[ibin],
                "acceptance_den_err": acc_den_err[ibin],
                "acceptance_ratio": ratio[ibin],
                "acceptance_ratio_err": ratio_err[ibin],
                "integrated_acceptance_num": integrated_num,
                "integrated_acceptance_num_err": integrated_num_err,
                "integrated_acceptance_den": integrated_den,
                "integrated_acceptance_den_err": integrated_den_err,
                "integrated_acceptance_ratio": integrated_ratio,
                "integrated_acceptance_ratio_err": integrated_ratio_err,
            })
        #endfor
    #endfor

    fig.suptitle(
        rf"ROOT-level acceptance ratio vs {cfg.title}",
        fontsize=16,
    )

    out_png = os.path.join(
        output_dir,
        f"acceptance_ratio_{cfg.key}.png",
    )

    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    return rows


def write_csv(path: str, rows: List[Dict]) -> None:
    if not rows:
        with open(path, "w") as f:
            f.write("")
        #endwith

        return
    #endif

    import csv

    fieldnames = list(rows[0].keys())

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)
        #endfor
    #endwith


def write_summary(path: str, file_results: Dict, acceptances: Dict, diagnostic_rows: List[Dict]) -> None:
    with open(path, "w") as f:
        f.write("ROOT-level DVCS MC acceptance comparison\n")
        f.write("========================================\n\n")

        f.write("Interpretation note\n")
        f.write("-------------------\n")
        f.write(
            "This script computes acceptance directly from ROOT trees as N_rec/N_gen.\n"
            "It does not apply the later DVCS analysis binning, topology filters,\n"
            "fiducial/exclusivity cuts, or CSV-level validity masks. Therefore, if\n"
            "these ROOT-level ratios are similar for inbending and outbending but the\n"
            "CSV-level ratios are not, the discrepancy is likely introduced downstream\n"
            "in the analysis code that reads and filters these trees.\n\n"
        )

        f.write("Files processed\n")
        f.write("---------------\n")

        for period in sorted(file_results.keys()):
            for sample_type in ["gen", "rec"]:
                r = file_results[period][sample_type]
                f.write(
                    f"{period:10s} {sample_type:3s} "
                    f"entries={r['n_entries_tree']} "
                    f"tree={r['tree']} "
                    f"path={r['path']}\n"
                )

                if r.get("missing_branches"):
                    f.write(f"  missing branches: {', '.join(r['missing_branches'])}\n")
                #endif
            #endfor
        #endfor

        f.write("\nTotal in-range counts by variable\n")
        f.write("---------------------------------\n")

        for cfg in VARIABLES:
            f.write(f"\n{cfg.key} ({cfg.branch})\n")

            for period in ["Fa18 Inb", "Sp18 Inb", "Fa18 Out", "Sp18 Out"]:
                if cfg.key not in acceptances[period]:
                    continue
                #endif

                gen_total = np.sum(acceptances[period][cfg.key]["gen_counts"])
                rec_total = np.sum(acceptances[period][cfg.key]["rec_counts"])
                acc_total = rec_total / gen_total if gen_total > 0.0 else np.nan

                f.write(
                    f"  {period:10s} gen={gen_total:12.0f} "
                    f"rec={rec_total:12.0f} rec/gen={acc_total:.6g}\n"
                )
            #endfor
        #endfor

        f.write("\nBranch finite ranges after transformations\n")
        f.write("------------------------------------------\n")

        for cfg in VARIABLES:
            f.write(f"\n{cfg.key} ({cfg.branch})\n")

            for period in ["Fa18 Inb", "Sp18 Inb", "Fa18 Out", "Sp18 Out"]:
                for sample_type in ["gen", "rec"]:
                    r = file_results[period][sample_type]

                    if cfg.key not in r["variables"]:
                        continue
                    #endif

                    vr = r["variables"][cfg.key]

                    f.write(
                        f"  {period:10s} {sample_type:3s} "
                        f"finite_min={vr['finite_min']:.6g} "
                        f"finite_max={vr['finite_max']:.6g} "
                        f"n_finite={vr['n_finite']} "
                        f"n_in_range={vr['n_in_range']}\n"
                    )
                #endfor
            #endfor
        #endfor

        f.write("\nMedian acceptance ratios by variable\n")
        f.write("------------------------------------\n")

        for cfg in VARIABLES:
            f.write(f"\n{cfg.key}\n")

            for ratio_label in [
                "Sp18 Inb / Fa18 Inb acceptance",
                "Sp18 Out / Fa18 Out acceptance",
            ]:
                vals = [
                    row["acceptance_ratio"]
                    for row in diagnostic_rows
                    if (
                        row["variable"] == cfg.key
                        and row["ratio_label"] == ratio_label
                        and np.isfinite(row["acceptance_ratio"])
                    )
                ]

                integrated_vals = [
                    row["integrated_acceptance_ratio"]
                    for row in diagnostic_rows
                    if (
                        row["variable"] == cfg.key
                        and row["ratio_label"] == ratio_label
                        and np.isfinite(row["integrated_acceptance_ratio"])
                    )
                ]

                if not vals:
                    continue
                #endif

                vals = np.asarray(vals, dtype=float)

                integrated_ratio = np.nan
                if integrated_vals:
                    integrated_ratio = float(integrated_vals[0])
                #endif

                f.write(
                    f"  {ratio_label:35s} "
                    f"median={np.nanmedian(vals):.4g} "
                    f"p16={np.nanpercentile(vals, 16):.4g} "
                    f"p84={np.nanpercentile(vals, 84):.4g} "
                    f"integrated={integrated_ratio:.4g}\n"
                )
            #endfor
        #endfor
    #endwith


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ROOT-level generated/reconstructed MC acceptance-ratio diagnostic."
    )
    parser.add_argument(
        "--outdir",
        default=OUTPUT_DIR_DEFAULT,
        help=f"Output directory. Default: {OUTPUT_DIR_DEFAULT}",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=MAX_WORKERS_DEFAULT,
        help=f"Maximum parallel workers. Default: {MAX_WORKERS_DEFAULT}",
    )
    parser.add_argument(
        "--tree",
        default=None,
        help="Optional TTree name. If omitted, first TTree in each file is used.",
    )
    parser.add_argument(
        "--step-size",
        default="200 MB",
        help="uproot iterate step size. Default: 200 MB.",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    tasks = []

    for period, paths in FILES.items():
        for sample_type in ["gen", "rec"]:
            tasks.append((period, sample_type, paths[sample_type]))
        #endfor
    #endfor

    file_results = {period: {} for period in FILES.keys()}

    print(f"[root acceptance] Processing {len(tasks)} files with max_workers={args.max_workers}")

    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        future_to_task = {}

        for period, sample_type, path in tasks:
            fut = executor.submit(
                hist_file_task,
                period,
                sample_type,
                path,
                args.tree,
                VARIABLES,
                args.step_size,
            )
            future_to_task[fut] = (period, sample_type, path)
        #endfor

        for fut in as_completed(future_to_task):
            period, sample_type, path = future_to_task[fut]

            try:
                result = fut.result()
            except Exception as exc:
                raise RuntimeError(
                    f"Failed while processing {period} {sample_type} {path}: {exc}"
                ) from exc
            #endtry

            file_results[period][sample_type] = result

            print(
                f"[root acceptance] Done {period:10s} {sample_type:3s}: "
                f"entries={result['n_entries_tree']} tree={result['tree']}"
            )
        #endfor
    #endwith

    acceptances = {period: {} for period in FILES.keys()}

    for period in FILES.keys():
        gen_vars = file_results[period]["gen"]["variables"]
        rec_vars = file_results[period]["rec"]["variables"]

        for cfg in VARIABLES:
            if cfg.key not in gen_vars or cfg.key not in rec_vars:
                continue
            #endif

            gen_counts = gen_vars[cfg.key]["counts"]
            rec_counts = rec_vars[cfg.key]["counts"]

            acceptance, acceptance_err = compute_acceptance(gen_counts, rec_counts)

            acceptances[period][cfg.key] = {
                "edges": gen_vars[cfg.key]["edges"],
                "gen_counts": gen_counts,
                "rec_counts": rec_counts,
                "acceptance": acceptance,
                "acceptance_err": acceptance_err,
            }
        #endfor
    #endfor

    all_rows = []

    for cfg in VARIABLES:
        available = all(
            cfg.key in acceptances[period]
            for period in ["Fa18 Inb", "Sp18 Inb", "Fa18 Out", "Sp18 Out"]
        )

        if not available:
            print(f"[root acceptance] Skipping {cfg.key}: missing in at least one file")
            continue
        #endif

        rows = make_variable_plot(cfg, acceptances, args.outdir)
        all_rows.extend(rows)

        print(f"[root acceptance] Wrote plot for {cfg.key}")
    #endfor

    diag_csv = os.path.join(args.outdir, "root_acceptance_ratio_diagnostics.csv")
    summary_txt = os.path.join(args.outdir, "root_acceptance_ratio_summary.txt")

    write_csv(diag_csv, all_rows)
    write_summary(summary_txt, file_results, acceptances, all_rows)

    print("\n[root acceptance] Wrote:")
    print(f"  {diag_csv}")
    print(f"  {summary_txt}")
    print(f"  PNGs in {args.outdir}")


if __name__ == "__main__":
    main()