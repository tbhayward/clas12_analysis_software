#!/usr/bin/env python3
import os
from load_binning_scheme import load_binning_scheme
from efficiency_analysis import (plot_mc_normalized_efficiencies, 
                                plot_data_normalized_efficiencies,
                                plot_DAF_data_normalized_efficiencies, 
                                 calculate_binned_efficiencies_all, 
                                 plot_binned_efficiencies_all)

def main():
    # --- Load binning scheme ---
    csv_file_path = os.path.join("imports", "integrated_bin_v2.csv")
    binning_scheme = load_binning_scheme(csv_file_path)
    print("Loaded binning scheme:")
    for b in binning_scheme:
        print(b)
    #endfor

    # --- (Optional) Exclusivity processing ---
    # The following is an example of how you might run your multi-stage exclusivity processing.
    # This code is commented out because you may have already run it or want to run it separately.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    periods_to_run = [
        ("DVCS_Fa18_inb",  "dvcs"),
        ("DVCS_Fa18_out",  "dvcs"),
        ("DVCS_Sp19_inb",  "dvcs"),
        ("eppi0_Fa18_inb", "eppi0"),
        ("eppi0_Fa18_out", "eppi0"),
        ("eppi0_Sp19_inb", "eppi0"),
    ]
    tasks = [(period, analysis_type, "exclusivity_output") for period, analysis_type in periods_to_run]
    with ProcessPoolExecutor(max_workers=6) as executor:
        futures = [executor.submit(run_period, task) for task in tasks]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                print(f"⚠️ A task generated an exception: {exc}")
            #endif
        #endfor
    print("🧩 Combining exclusivity results (JSON files from each topology & stage)...")
    combine_results("exclusivity_output")
    """
    #endif

    # --- Analyze and plot overall normalized efficiency ---
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # plot_mc_normalized_efficiencies(output_dir)
    # plot_data_normalized_efficiencies(output_dir)
    plot_DAF_data_normalized_efficiencies(output_dir)
    
    # # --- Calculate binned efficiencies for all run periods ---
    # # This function loops over all run periods and calculates binned efficiencies,
    # # saving one JSON file per run period.
    # calculate_binned_efficiencies_all(binning_scheme, output_dir)
    
    # # --- Plot binned efficiencies for all run periods ---
    # plot_binned_efficiencies_all(binning_csv=csv_file_path,
    #                              output_dir=os.path.join(output_dir, "efficiency_plots"))
#enddef

if __name__ == "__main__":
    main()
#enddef