#!/usr/bin/env python3

import os
from concurrent.futures import ProcessPoolExecutor, as_completed

# We import the multi-stage process_period function and the combine_results function
from label_formatting import configure_global_style
from exclusivity import process_period_multi_stage, combine_results
from load_binning_scheme import load_binning_scheme

def run_period(args):
    """
    Helper function for parallel processing.
    Each process calls process_period_multi_stage for one (period, analysis_type) tuple.
    """
    period, analysis_type, output_dir = args
    # It is a good idea to configure ROOT style separately in each process if needed.
    configure_global_style()
    process_period_multi_stage(period, output_dir, analysis_type)

def main():
    """
    The 'main' entry point.
    - Applies global style (optional).
    - Creates output directory.
    - Loops over run periods (DVCS or eppi0).
    - Calls process_period_multi_stage(...) for each in parallel (3 at a time).
    - Combines final JSON outputs if you like.
    """
    # Optional: configure global ROOT style once at startup
    configure_global_style()
    
    output_dir = "exclusivity"
    os.makedirs(output_dir, exist_ok=True)

    print("🚀 Starting multi-stage analysis with parallel processing\n")

    # EXAMPLE: Suppose you want to process DVCS for Fa18_inb, Fa18_out, Sp19_inb,
    # and also eppi0 for some set of periods.
    periods_to_run = [
        ("DVCS_Fa18_inb",  "dvcs"),
        ("DVCS_Fa18_out",  "dvcs"),
        ("DVCS_Sp19_inb",  "dvcs"),
        ("eppi0_Fa18_inb", "eppi0"),
        ("eppi0_Fa18_out", "eppi0"),
        ("eppi0_Sp19_inb", "eppi0"),
    ]

    # # Prepare a list of argument tuples for each period (include the output_dir)
    # tasks = [(period, analysis_type, output_dir) for period, analysis_type in periods_to_run]

    # # Use ProcessPoolExecutor with a maximum of 3 workers
    # with ProcessPoolExecutor(max_workers=6) as executor:
    #     futures = [executor.submit(run_period, task) for task in tasks]
    #     for future in as_completed(futures):
    #         try:
    #             future.result()  # wait for individual task to complete
    #         except Exception as exc:
    #             print(f"⚠️ A task generated an exception: {exc}")
    # print("🧩 Combining results (JSON files from each topology & stage)...")
    # combine_results(output_dir)


    # Load the binning scheme from the CSV file.
    csv_file_path = os.path.join("imports", "integrated_bin_v2.csv")
    binning_scheme = load_binning_scheme(csv_file_path)
    print("Loaded binning scheme:")
    for b in binning_scheme:
        print(b)
    print("\n🎉 Analysis complete!")

if __name__ == "__main__":
    main()