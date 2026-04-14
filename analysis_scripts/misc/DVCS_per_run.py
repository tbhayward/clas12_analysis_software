import uproot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import defaultdict

# Create output directory if it doesn't exist
output_dir = "output/dvcs_per_run"
os.makedirs(output_dir, exist_ok=True)

# Path to CSV with run information
csv_file = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity/global.csv"

# Read CSV, skipping comment lines
run_info_df = pd.read_csv(
    csv_file,
    comment="#",
    header=None,
    names=["runnum", "charge_nC", "col2", "col3", "col4", "col5"],
)

# Create a dictionary mapping run number to accumulated charge
run_charge_map = dict(zip(run_info_df["runnum"], run_info_df["charge_nC"]))

# -----------------------
# Run -> current maps
# (ported from yield_totals.cpp)
# -----------------------

FA18_INB_CURRENT = {
    # 5 nA
    5418: 5, 5419: 5,

    # 40 nA
    5335: 40, 5339: 40, 5341: 40,
    5340: 40, 5342: 40, 5343: 40, 5344: 40,

    # 45 nA
    5032: 45, 5036: 45, 5038: 45, 5039: 45, 5040: 45, 5041: 45,
    5043: 45, 5045: 45, 5052: 45, 5053: 45, 5116: 45, 5117: 45,
    5119: 45, 5120: 45, 5124: 45, 5125: 45, 5126: 45, 5127: 45,
    5139: 45, 5153: 45, 5158: 45, 5162: 45, 5163: 45, 5164: 45,
    5181: 45, 5191: 45, 5193: 45, 5195: 45, 5196: 45, 5197: 45,
    5198: 45, 5199: 45, 5200: 45, 5201: 45, 5202: 45, 5203: 45,
    5204: 45, 5205: 45, 5206: 45, 5208: 45, 5211: 45, 5212: 45,
    5215: 45, 5216: 45, 5219: 45, 5220: 45, 5221: 45, 5222: 45,
    5223: 45, 5230: 45, 5231: 45, 5232: 45, 5233: 45, 5234: 45,
    5235: 45, 5237: 45, 5238: 45, 5247: 45, 5248: 45, 5249: 45,
    5252: 45, 5253: 45, 5257: 45, 5258: 45, 5259: 45, 5261: 45,
    5262: 45, 5303: 45, 5304: 45, 5305: 45, 5306: 45, 5307: 45,
    5310: 45, 5311: 45, 5315: 45, 5317: 45, 5318: 45, 5319: 45,
    5320: 45, 5323: 45, 5324: 45, 5333: 45, 5334: 45, 5346: 45,
    5347: 45, 5349: 45, 5351: 45, 5354: 45, 5355: 45, 5367: 45,

    # Extra Fa18 Inb 45 nA runs
    5046: 45, 5047: 45, 5051: 45,
    5128: 45, 5129: 45, 5130: 45,
    5159: 45, 5160: 45,
    5165: 45, 5166: 45, 5167: 45, 5168: 45, 5169: 45,
    5180: 45, 5182: 45, 5183: 45,
    5190: 45,
    5239: 45,
    5336: 45,

    # 50 nA
    5356: 50, 5357: 50, 5358: 50, 5359: 50, 5360: 50, 5361: 50,
    5362: 50, 5366: 50,

    # 55 nA
    5368: 55, 5369: 55, 5372: 55, 5373: 55, 5374: 55, 5375: 55,
    5376: 55, 5377: 55, 5378: 55, 5379: 55, 5380: 55, 5381: 55,
    5382: 55, 5383: 55, 5386: 55, 5390: 55, 5391: 55, 5392: 55,
    5393: 55, 5398: 55, 5400: 55, 5401: 55, 5403: 55, 5404: 55,
    5406: 55, 5407: 55,
}

FA18_OUT_CURRENT = {
    # 5 nA
    5443: 5,

    # 20 nA
    5444: 20,

    # 40 nA
    5423: 40, 5424: 40, 5425: 40, 5426: 40, 5428: 40, 5429: 40,
    5430: 40, 5432: 40, 5434: 40, 5435: 40, 5436: 40, 5437: 40,
    5438: 40, 5440: 40, 5441: 40, 5442: 40, 5445: 40, 5447: 40,
    5449: 40, 5450: 40, 5451: 40, 5452: 40, 5453: 40, 5454: 40,
    5455: 40, 5460: 40, 5464: 40, 5465: 40, 5466: 40, 5467: 40,
    5468: 40, 5469: 40, 5470: 40, 5471: 40, 5472: 40, 5473: 40,
    5474: 40, 5475: 40, 5476: 40, 5478: 40, 5479: 40, 5480: 40,
    5481: 40, 5482: 40, 5483: 40, 5485: 40, 5486: 40, 5487: 40,
    5497: 40, 5498: 40, 5499: 40, 5500: 40, 5504: 40,

    # Extra Fa18 Out 40 nA runs
    5448: 40, 5495: 40, 5496: 40,

    # 50 nA
    5507: 50, 5516: 50, 5517: 50, 5518: 50, 5519: 50, 5520: 50,
    5521: 50, 5522: 50, 5523: 50, 5524: 50, 5525: 50, 5526: 50,
    5527: 50, 5528: 50, 5530: 50, 5532: 50, 5533: 50, 5534: 50,
    5535: 50, 5536: 50, 5537: 50, 5538: 50, 5540: 50, 5541: 50,
    5543: 50, 5544: 50, 5545: 50, 5546: 50, 5547: 50, 5548: 50,
    5549: 50, 5550: 50, 5551: 50, 5552: 50, 5555: 50, 5556: 50,
    5557: 50, 5558: 50, 5559: 50, 5562: 50, 5569: 50, 5570: 50,
    5571: 50, 5572: 50, 5573: 50, 5574: 50, 5577: 50, 5578: 50,
    5591: 50, 5592: 50, 5594: 50, 5597: 50, 5598: 50, 5600: 50,
    5601: 50, 5602: 50, 5603: 50, 5604: 50, 5606: 50, 5607: 50,
    5611: 50, 5612: 50, 5613: 50, 5614: 50, 5616: 50, 5618: 50,
    5619: 50, 5624: 50, 5625: 50, 5626: 50, 5627: 50, 5628: 50,
    5629: 50, 5630: 50, 5631: 50, 5632: 50, 5633: 50, 5635: 50,
    5637: 50, 5638: 50, 5639: 50, 5641: 50, 5643: 50, 5644: 50,
    5645: 50, 5646: 50, 5647: 50, 5648: 50, 5649: 50, 5650: 50,
    5651: 50, 5652: 50, 5654: 50, 5655: 50, 5656: 50, 5662: 50,
    5663: 50, 5664: 50, 5665: 50, 5666: 50,
    5615: 50,

    # Extra Fa18 Out 50 nA runs
    5505: 50, 5567: 50, 5617: 50, 5621: 50, 5623: 50,

    # Extra Fa18 Out 5 nA run
    5610: 5,
}


def resolve_current(period_label, runnum):
    """
    Resolve beam current (nA) for a given period and run number.

    Returns (True, current_nA) if known, else (False, None).
    """
    label = period_label.lower()

    # Fa18 Inb
    if label == "rga_fa18_inb":
        if runnum in FA18_INB_CURRENT:
            return True, FA18_INB_CURRENT[runnum]
        #endif
        return False, None
    #endif

    # Fa18 Out
    if label == "rga_fa18_out":
        if runnum in FA18_OUT_CURRENT:
            return True, FA18_OUT_CURRENT[runnum]
        #endif
        return False, None
    #endif

    # Sp18 Out: 30 nA (3211-3293), 45 nA (3867-3987)
    if label == "rga_sp18_out":
        if 3211 <= runnum <= 3293:
            return True, 30
        #endif
        if 3867 <= runnum <= 3987:
            return True, 45
        #endif
        return False, None
    #endif

    # Sp18 Inb: overrides, then ranges
    if label == "rga_sp18_inb":
        if runnum == 3418:
            return True, 70
        #endif
        if runnum == 3421 or runnum == 3422:
            return True, 35
        #endif
        if runnum == 3429:
            return True, 50
        #endif

        # 35 nA from 3306-3411
        if 3306 <= runnum <= 3411:
            return True, 35
        #endif
        # 50 nA from 3431-4325
        if 3431 <= runnum <= 4325:
            return True, 50
        #endif

        return False, None
    #endif

    # Sp19 Inb: all 50 nA
    if label == "rga_sp19_inb":
        if runnum == 6616:
            return True, 5
        #endif
        return True, 50
    #endif

    # Other periods: unknown
    return False, None
#enddef


# Define all periods and their ROOT files, with fa18_inb first
period_files = [
    ("rga_fa18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root"),
    ("rga_fa18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root"),
    ("rga_sp19_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"),
    ("rga_sp18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root"),
    ("rga_sp18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root"),
]


for period_label, root_path in period_files:
    print("\n" + "=" * 80)
    print(f"Processing period: {period_label}")
    print("=" * 80)

    # Check that the ROOT file exists
    if not os.path.exists(root_path):
        print(f"Warning: ROOT file for {period_label} not found at {root_path}, skipping.")
        #endif
        continue
    #endif

    # Open the ROOT file and access the PhysicsEvents tree
    root_file = uproot.open(root_path)
    if "PhysicsEvents" not in root_file:
        print(f"Warning: 'PhysicsEvents' tree not found in {root_path}, skipping.")
        #endif
        continue
    #endif

    tree = root_file["PhysicsEvents"]

    # Get the runnum branch data
    runnum_data = tree["runnum"].array(library="np")

    # Count events per run number
    unique_runs, event_counts = np.unique(runnum_data, return_counts=True)

    # Group by current: current_nA -> lists of runs / values / errors
    per_current_runs = defaultdict(list)
    per_current_values = defaultdict(list)
    per_current_errors = defaultdict(list)

    missing_charge_count = 0
    zero_charge_count = 0
    unknown_current_count = 0

    # Loop over runs in this period
    for run_num, count in zip(unique_runs, event_counts):
        run_num_int = int(run_num)
        count_int = int(count)

        if run_num_int not in run_charge_map:
            missing_charge_count += 1
            continue
        #endif

        charge = run_charge_map[run_num_int]
        if charge <= 0.0:
            print(f"Warning: Run {run_num_int} has zero or negative charge, skipping.")
            zero_charge_count += 1
            continue
        #endif

        # Events per nC and stat uncertainty
        value = count_int / charge
        error = np.sqrt(count_int) / charge if count_int > 0.0 else 0.0

        # Resolve current
        ok_cur, current_nA = resolve_current(period_label, run_num_int)
        if not ok_cur:
            print(f"Warning: Run {run_num_int} has charge info but no current mapping for period {period_label}; skipping in plot.")
            unknown_current_count += 1
            continue
        #endif

        per_current_runs[current_nA].append(run_num_int)
        per_current_values[current_nA].append(value)
        per_current_errors[current_nA].append(error)
    #endfor

    # If no runs with both charge and current info, skip plotting
    if not per_current_runs:
        print(f"\nNo runs with both charge and current mapping were found for period {period_label}.")
        if missing_charge_count > 0:
            print(f"  Note: {missing_charge_count} runs from this ROOT file are missing in the charge CSV.")
        #endif
        if zero_charge_count > 0:
            print(f"  Note: {zero_charge_count} runs had zero or negative charge recorded.")
        #endif
        if unknown_current_count > 0:
            print(f"  Note: {unknown_current_count} runs had no current mapping and were skipped.")
        #endif
        print("Skipping plot for this period.\n")
        #endif
        continue
    #endif

    # Compute statistics per current and prepare for plotting
    per_current_stats = {}
    total_valid_runs = 0

    for current, run_list in per_current_runs.items():
        runs_arr = np.array(run_list)
        vals_arr = np.array(per_current_values[current])
        errs_arr = np.array(per_current_errors[current])

        # Sort by run number
        sort_idx = np.argsort(runs_arr)
        runs_sorted = runs_arr[sort_idx]
        vals_sorted = vals_arr[sort_idx]
        errs_sorted = errs_arr[sort_idx]

        # -------------------------
        # Stage 1: initial mean and >10 sigma (per-run) extreme outliers
        # -------------------------
        if len(vals_sorted) > 1:
            mean0 = np.mean(vals_sorted)
        else:
            mean0 = vals_sorted[0]
        #endif

        extreme_mask = np.abs(vals_sorted - mean0) > 10.0 * errs_sorted

        # If all points are extreme, fall back to no trimming
        if np.any(~extreme_mask):
            trimmed_vals = vals_sorted[~extreme_mask]
        else:
            trimmed_vals = vals_sorted
            extreme_mask = np.zeros_like(extreme_mask, dtype=bool)
        #endif

        # -------------------------
        # Stage 2: trimmed mean and 5 sigma (per-run) outliers
        # -------------------------
        if len(trimmed_vals) > 1:
            mean_trimmed = np.mean(trimmed_vals)
            sigma_spread = np.std(trimmed_vals, ddof=1)
        else:
            mean_trimmed = trimmed_vals[0]
            sigma_spread = 0.0
        #endif

        mild_mask = np.abs(vals_sorted - mean_trimmed) > 5.0 * errs_sorted
        final_outlier_mask = extreme_mask | mild_mask

        per_current_stats[current] = {
            "runs": runs_sorted,
            "vals": vals_sorted,
            "errs": errs_sorted,
            "mean": mean_trimmed,
            "sigma_spread": sigma_spread,
            "extreme_mask": extreme_mask,
            "outliers": final_outlier_mask,
        }
        total_valid_runs += len(runs_sorted)
    #endfor

    # Print statistics
    print(f"\nStatistics for period {period_label}:")
    for current in sorted(per_current_stats.keys()):
        st = per_current_stats[current]
        print(
            f"  Current {current} nA: "
            f"mean_trimmed(events/nC) = {st['mean']:.6f}, "
            f"sigma(spread, trimmed) = {st['sigma_spread']:.6f}, "
            f"N_runs = {len(st['runs'])}"
        )
    #endfor

    # Extra diagnostic for Sp18 Out 30 nA
    if period_label == "rga_sp18_out" and 30 in per_current_stats:
        print("\nFULL 30 nA run list for rga_sp18_out:")
        st = per_current_stats[30]
        for run, val, err in zip(st["runs"], st["vals"], st["errs"]):
            print(f"  run {int(run)}: {val:.6f} +/- {err:.6f}")
        #endfor

        weighted_num = np.sum(st["vals"] / (st["errs"] * st["errs"]))
        weighted_den = np.sum(1.0 / (st["errs"] * st["errs"]))
        weighted_mean = weighted_num / weighted_den if weighted_den > 0.0 else float("nan")

        print("")
        print(f"  Unweighted trimmed mean shown by dashed line: {st['mean']:.6f}")
        print(f"  Inverse-variance weighted mean of visible 30 nA runs: {weighted_mean:.6f}")

        total_counts_30 = 0.0
        total_charge_30 = 0.0
        for run, val in zip(st["runs"], st["vals"]):
            run_int = int(run)
            if run_int in run_charge_map:
                charge = float(run_charge_map[run_int])
                counts = val * charge
                total_counts_30 += counts
                total_charge_30 += charge
            #endif
        #endfor

        aggregate_rate_30 = total_counts_30 / total_charge_30 if total_charge_30 > 0.0 else float("nan")
        print(f"  Aggregate sum(counts)/sum(charge) for 30 nA runs in this script: {aggregate_rate_30:.6f}")
    #endif

    # Print outliers
    print("\nRuns more than 5 sigma (per-run stat error) from the trimmed mean (and >10 sigma extremes):")
    any_outliers = False
    for current in sorted(per_current_stats.keys()):
        st = per_current_stats[current]
        runs = st["runs"]
        vals = st["vals"]
        errs = st["errs"]
        outmask = st["outliers"]
        extreme_mask = st["extreme_mask"]

        for run, val, err, is_out, is_extreme in zip(runs, vals, errs, outmask, extreme_mask):
            if is_out:
                if not any_outliers:
                    any_outliers = True
                #endif
                tag = "EXTREME(>10sigma)" if is_extreme else ">5sigma"
                print(f"  Run {int(run)} (current {current} nA): {val:.6f} +/- {err:.6f} events/nC  [{tag}]")
            #endif
        #endfor
    #endfor
    if not any_outliers:
        print("  None")
    #endif

    # Create plot
    plt.figure(figsize=(12, 6))

    # Color map per current (sorted for consistency)
    color_cycle = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
    sorted_currents = sorted(per_current_stats.keys())
    current_colors = {}
    for idx, current in enumerate(sorted_currents):
        current_colors[current] = color_cycle[idx % len(color_cycle)]
    #endfor

    # Plot per current: inliers as solid circles, outliers as open circles
    for current in sorted_currents:
        st = per_current_stats[current]
        runs = st["runs"]
        vals = st["vals"]
        errs = st["errs"]
        outmask = st["outliers"]
        color = current_colors[current]

        inmask = ~outmask

        # Inliers (solid circles) with label for legend
        if np.any(inmask):
            plt.errorbar(
                runs[inmask],
                vals[inmask],
                yerr=errs[inmask],
                fmt="o",
                markersize=4,
                linestyle="none",
                markerfacecolor=color,
                markeredgecolor=color,
                ecolor=color,
                label=f"{current} nA",
            )
        #endif

        # Outliers (open circles), no label (same color)
        if np.any(outmask):
            plt.errorbar(
                runs[outmask],
                vals[outmask],
                yerr=errs[outmask],
                fmt="o",
                markersize=6,
                linestyle="none",
                markerfacecolor="none",
                markeredgecolor=color,
                ecolor=color,
            )
        #endif

        # Mean line for this current (trimmed mean)
        plt.axhline(st["mean"], color=color, linestyle="--", linewidth=1.0)
    #endfor

    plt.xlabel("Run Number")
    plt.ylabel("Events / nC")
    plt.title(f"DVCS Events per nC by Run Number ({period_label})")
    plt.grid(True, alpha=0.3)

    # Force the same y-axis range on all plots
    plt.ylim(0.0, 0.10)

    plt.xticks(rotation=45)
    plt.legend(title="Beam current (nA)", fontsize=10)
    plt.tight_layout()

    # Save the plot
    out_path = os.path.join(output_dir, f"dvcs_per_nC_{period_label}.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"\nPlot saved to {out_path}")
    print(f"Processed {total_valid_runs} runs with valid charge and current info for {period_label}")
    if missing_charge_count > 0:
        print(f"{missing_charge_count} runs were missing from the charge CSV for this period.")
    #endif
    if zero_charge_count > 0:
        print(f"{zero_charge_count} runs had zero or negative charge recorded and were skipped.")
    #endif
    if unknown_current_count > 0:
        print(f"{unknown_current_count} runs had no current mapping and were skipped.")
    #endif
#endfor