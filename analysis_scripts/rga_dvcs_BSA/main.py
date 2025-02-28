#!/usr/bin/env python3

import os
import json
from concurrent.futures import ProcessPoolExecutor, as_completed

# We import the multi-stage process_period function and combine_results from exclusivity.
from label_formatting import configure_global_style
from exclusivity import process_period_multi_stage, combine_results
from load_binning_scheme import load_binning_scheme
from calculate_contamination import calculate_contamination

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
    # print("Loaded binning scheme:")
    # for b in binning_scheme:
    #     print(b)
    print(f"Loaded binning_scheme: {type(binning_scheme)} - {binning_scheme}")
    
    # --- Contamination calculation ---
    # Define DVCS periods for contamination calculation.
    dvcs_periods = [
        ("DVCS_Fa18_inb", "dvcs"),
        ("DVCS_Fa18_out", "dvcs"),
        ("DVCS_Sp19_inb", "dvcs")
    ]
    # Define the three topologies.
    topologies = ["(FD,FD)", "(CD,FD)", "(CD,FT)"]
    
    # Build tasks: each task is (period, topology, analysis_type, binning_scheme)
    tasks = []
    for period, analysis_type in dvcs_periods:
        for topo in topologies:
            tasks.append((period, topo, analysis_type, binning_scheme))
    
    # print("HELLO WORLD 1")
    # # Run contamination calculations in parallel (max 3 workers).
    # with ProcessPoolExecutor(max_workers=1) as executor:
    #     future_to_task = {executor.submit(calculate_contamination, *task): task for task in tasks}
    #     for future in as_completed(future_to_task):
    #         print("HELLO WORLD 2")
    #         task = future_to_task[future]
    #         try:
    #             print("HELLO WORLD 3")
    #             period, topology, analysis_type, _ = task  # Unpack all 4 elements.
    #             print(period)
    #             print(topology)
    #             print(analysis_type)
    #             result = future.result()
    #             safe_topo = topology.replace("(", "").replace(")", "").replace(",", "_")
    #             print("HELLO WORLD 4")
    #             json_filename = f"contamination_{period}_{safe_topo}.json"
    #             json_path = os.path.join(contamination_dir, json_filename)
    #             with open(json_path, "w") as f:
    #                 json.dump(result, f, indent=2)
    #             print(f"Saved contamination for {period} {topology} to {json_path}")
    #         except Exception as exc:
    #             print(f"Task {task} generated an exception: {exc}")
    
    # print("\n🎉 Analysis complete!")

    print("HELLO WORLD 1")

    # Sequential execution for debugging
    for task in tasks:
        print("HELLO WORLD 2")
        try:
            print("HELLO WORLD 3")
            period, topology, analysis_type, _ = task  # Unpack all 4 elements.
            print(period)
            print(topology)
            print(analysis_type)
            
            print(f"Calling calculate_contamination with args: {task}")
            # Directly call `calculate_contamination` instead of using executor.submit
            result = calculate_contamination(*task)

            safe_topo = topology.replace("(", "").replace(")", "").replace(",", "_")
            print("HELLO WORLD 4")

            json_filename = f"contamination_{period}_{safe_topo}.json"
            json_path = os.path.join(contamination_dir, json_filename)

            with open(json_path, "w") as f:
                json.dump(result, f, indent=2)

            print(f"Saved contamination for {period} {topology} to {json_path}")

        except Exception as exc:
            print(f"Task {task} generated an exception: {exc}")

    print("\n🎉 Analysis complete!")

if __name__ == "__main__":
    main()