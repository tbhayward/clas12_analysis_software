#!/usr/bin/env python3

import os
import json
import math
from concurrent.futures import ProcessPoolExecutor, as_completed

# We import the multi-stage process_period function and combine_results from exclusivity.
from label_formatting import configure_global_style
from exclusivity import process_period_multi_stage, combine_results
from load_binning_scheme import load_binning_scheme
from calculate_contamination import calculate_contamination
from calculate_bin_means import calculate_bin_means
from plot_contamination import plot_contamination_for_run
from determine_raw_bsa import calculate_raw_bsa
from determine_final_bsa import determine_final_bsa
from plot_bsa import plot_raw_bsa, plot_adjusted_bsa, plot_combined_bsa

def run_period(args):
    """
    Helper function for parallel processing.
    Each process calls process_period_multi_stage for one (period, analysis_type) tuple.
    """
    period, analysis_type, output_dir = args
    configure_global_style()
    process_period_multi_stage(period, output_dir, analysis_type)

def main():
    """
    Main entry point.
    - Configures global ROOT style.
    - Creates output directories.
    - Optionally runs the multi-stage exclusivity processing in parallel (currently commented out).
    - Loads the binning scheme from the CSV file.
    - For each DVCS period and topology, calls calculate_contamination to compute 4D contamination,
      using the cuts from the combined cuts JSON.
    - Saves each contamination result as a JSON file in the "contamination" directory.
    """
    configure_global_style()
    
    output_dir = "exclusivity"
    os.makedirs(output_dir, exist_ok=True)
    contamination_dir = "contamination"
    os.makedirs(contamination_dir, exist_ok=True)
    
    print("🚀 Starting multi-stage analysis with parallel processing\n")
    
    # # --- Exclusivity processing (commented out) ---
    # periods_to_run = [
    #     ("DVCS_Fa18_inb",  "dvcs"),
    #     ("DVCS_Fa18_out",  "dvcs"),
    #     ("DVCS_Sp19_inb",  "dvcs"),
    #     ("eppi0_Fa18_inb", "eppi0"),
    #     ("eppi0_Fa18_out", "eppi0"),
    #     ("eppi0_Sp19_inb", "eppi0"),
    # ]
    # tasks = [(period, analysis_type, output_dir) for period, analysis_type in periods_to_run]
    # with ProcessPoolExecutor(max_workers=6) as executor:
    #     futures = [executor.submit(run_period, task) for task in tasks]
    #     for future in as_completed(futures):
    #         try:
    #             future.result()
    #         except Exception as exc:
    #             print(f"⚠️ A task generated an exception: {exc}")
    # print("🧩 Combining exclusivity results (JSON files from each topology & stage)...")
    # combine_results(output_dir)
    
    # --- Load binning scheme ---
    csv_file_path = os.path.join("imports", "integrated_bin_v2.csv")
    binning_scheme = load_binning_scheme(csv_file_path)
    print("Loaded binning scheme:")
    for b in binning_scheme:
        print(b)

    # calculate global means of bins
    # Define the DVCS periods (which we want to combine)
    dvcs_periods = [
        "DVCS_Fa18_inb",
        "DVCS_Fa18_out",
        "DVCS_Sp19_inb"
    ]
    # Define topologies to include
    topologies = ["(FD,FD)", "(CD,FD)", "(CD,FT)"]
    #endfor
    # Set your analysis type (usually "dvcs") 
    analysis_type = "dvcs"
    # Define where you want the final JSON to be saved
    output_json = "bin_means_global.json"
    # Call the updated function to calculate the GLOBAL bin means
    calculate_bin_means(dvcs_periods, topologies, analysis_type, binning_scheme, output_json)
    print("[main] Global bin-averaged kinematics have been computed!")

    
    # # --- Contamination calculation tasks ---
    # dvcs_periods = [
    #     ("DVCS_Fa18_inb", "dvcs"),
    #     ("DVCS_Fa18_out", "dvcs"),
    #     ("DVCS_Sp19_inb", "dvcs")
    # ]
    # topologies = ["(FD,FD)", "(CD,FD)", "(CD,FT)"]
    
    # # Build tasks: each task is (period, topology, analysis_type, binning_scheme)
    # tasks = []
    # for period, analysis_type in dvcs_periods:
    #     for topo in topologies:
    #         tasks.append((period, topo, analysis_type, binning_scheme))
    
    # # Run contamination calculations in parallel
    # period_results = {}  # Will group results by period
    # with ProcessPoolExecutor(max_workers=8) as executor:
    #     future_to_task = {executor.submit(calculate_contamination, *task): task for task in tasks}
    #     for future in as_completed(future_to_task):
    #         task = future_to_task[future]
    #         try:
    #             period, topology, analysis_type, _ = task
    #             result = future.result()
    #             # Group the result under the corresponding period:
    #             if period not in period_results:
    #                 period_results[period] = result
    #             else:
    #                 # Sum the counts in each bin (the keys should match if you use unique bins)
    #                 for key, counts in result.items():
    #                     if key not in period_results[period]:
    #                         period_results[period][key] = counts
    #                     else:
    #                         period_results[period][key]['N_data'] += counts['N_data']
    #                         period_results[period][key]['N_pi0_mc'] += counts['N_pi0_mc']
    #                         period_results[period][key]['N_pi0_exp'] += counts['N_pi0_exp']
    #                         period_results[period][key]['N_pi0_reco'] += counts['N_pi0_reco']
    #         except Exception as exc:
    #             print(f"Task {task} generated an exception: {exc}")
    
    # # Now recompute the contamination values for each bin from the summed counts,
    # # and write one file per period.
    # for period, result in period_results.items():
    #     for key, counts in result.items():
    #         N_data = counts['N_data']
    #         if N_data == 0 or counts['N_pi0_reco'] == 0:
    #             counts['c_i'] = 0.0
    #             counts['c_i_err'] = 0.0
    #         else:
    #             ratio = counts['N_pi0_exp'] / counts['N_pi0_reco']
    #             c_i = counts['N_pi0_mc'] * ratio / N_data
    #             rel_pi0_mc = 1 / math.sqrt(counts['N_pi0_mc']) if counts['N_pi0_mc'] > 0 else 0
    #             rel_pi0_exp = 1 / math.sqrt(counts['N_pi0_exp']) if counts['N_pi0_exp'] > 0 else 0
    #             rel_pi0_reco = 1 / math.sqrt(counts['N_pi0_reco']) if counts['N_pi0_reco'] > 0 else 0
    #             rel_data = 1 / math.sqrt(N_data)
    #             # Here we add the uncertainties for the ratio in quadrature.
    #             rel_ratio = math.sqrt(rel_pi0_exp**2 + rel_pi0_reco**2)
    #             rel_err = math.sqrt(rel_pi0_mc**2 + rel_ratio**2 + rel_data**2)
    #             c_i_err = c_i * rel_err
    #             counts['c_i'] = c_i
    #             counts['c_i_err'] = c_i_err
        
    #     # Filter out bins with zero contamination and round the numbers
    #     filtered_results = {
    #         str(key): {
    #             'c_i': round(counts['c_i'], 5),
    #             'c_i_err': round(counts['c_i_err'], 5)
    #         }
    #         for key, counts in result.items() if counts['c_i'] != 0
    #     }
        
    #     json_filename = f"contamination_{period}.json"
    #     json_path = os.path.join(contamination_dir, json_filename)
    #     with open(json_path, "w") as f:
    #         json.dump(filtered_results, f, indent=2)
    #     print(f"Saved combined contamination for {period} to {json_path}")

    # # --- Plotting contamination for each run period ---
    # run_periods = ["DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb"]
    # plots_dir = os.path.join("contamination", "contamination_plots")
    # os.makedirs(plots_dir, exist_ok=True)
    # csv_file_path = os.path.join("imports", "integrated_bin_v2.csv")

    # with ProcessPoolExecutor(max_workers=8) as executor:
    #     future_to_rp = {executor.submit(plot_contamination_for_run,
    #                                     run_period=rp,
    #                                     binning_csv=csv_file_path,
    #                                     contamination_dir="contamination",
    #                                     output_dir=plots_dir): rp for rp in run_periods}
    #     for future in as_completed(future_to_rp):
    #         rp = future_to_rp[future]
    #         try:
    #             future.result()
    #             print(f"Finished plotting contamination for {rp}")
    #         except Exception as exc:
    #             print(f"Plotting for {rp} generated an exception: {exc}")

        
    # # --- Raw BSA calculation ---
    # print("\n🚀 Calculating raw BSA values...")
    # csv_path = os.path.join("imports", "integrated_bin_v2.csv")
    # bsa_output = os.path.join("bsa_results")

    # # Create tasks using EXACT file_map keys
    # bsa_tasks = [
    #     # DVCS periods
    #     ("DVCS_Fa18_inb", "dvcs", csv_path, bsa_output),
    #     ("DVCS_Fa18_out", "dvcs", csv_path, bsa_output),
    #     ("DVCS_Sp19_inb", "dvcs", csv_path, bsa_output),
    #     # eppi0 periods
    #     ("eppi0_Fa18_inb", "eppi0", csv_path, bsa_output),
    #     ("eppi0_Fa18_out", "eppi0", csv_path, bsa_output),
    #     ("eppi0_Sp19_inb", "eppi0", csv_path, bsa_output)
    # ]

    # with ProcessPoolExecutor(max_workers=6) as executor:
    #     futures = {executor.submit(calculate_raw_bsa, *task): task for task in bsa_tasks}
    #     for future in as_completed(futures):
    #         task = futures[future]
    #         try:
    #             future.result()
    #             print(f"Finished BSA for {task[0]}")
    #         except Exception as exc:
    #             print(f"BSA failed for {task[0]}: {exc}")

    # # --- Final BSA calculation and combination ---
    # print("\n🔧 Calculating final adjusted BSA values...")
    # determine_final_bsa(
    #     contamination_dir="contamination",
    #     bsa_dir="bsa_results",
    #     final_dir="final_results"
    # )
    # print("✅ Final BSA results saved to final_results/ directory")

     # --- Plotting ---
    print("\n📊 Generating BSA plots...")
    csv_path = os.path.join("imports", "integrated_bin_v2.csv")
    
    plot_raw_bsa(csv_path)
    plot_adjusted_bsa(csv_path)
    plot_combined_bsa(csv_path)
    
    print("✅ All plots saved to bsa_plots/ directory")

    print("\n🎉 Analysis complete!")

if __name__ == "__main__":
    main()