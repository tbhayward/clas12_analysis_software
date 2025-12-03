import uproot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Create output directory if it doesn't exist
os.makedirs('output', exist_ok=True)

# Path to CSV with run information
csv_file = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Read CSV, skipping comment lines
run_info_df = pd.read_csv(
    csv_file,
    comment='#',
    header=None,
    names=['runnum', 'charge_nC', 'col2', 'col3', 'col4', 'col5']
)

# Create a dictionary mapping run number to accumulated charge
run_charge_map = dict(zip(run_info_df['runnum'], run_info_df['charge_nC']))

# Define all periods and their ROOT files
period_files = {
    "rga_fa18_out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
    "rga_fa18_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
    "rga_sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root",
    "rga_sp18_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root",
    "rga_sp18_out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root",
}

for period_label, root_path in period_files.items():
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
    runnum_data = tree["runnum"].array()

    # Count events per run number
    unique_runs, event_counts = np.unique(runnum_data, return_counts=True)

    # Calculate events per nC for each run
    events_per_nC = []
    valid_runs = []
    missing_charge_count = 0
    zero_charge_count = 0

    print("Runs with events/nC < 0.08:")
    print("-" * 30)

    for run_num, count in zip(unique_runs, event_counts):
        if run_num in run_charge_map:
            charge = run_charge_map[run_num]
            if charge > 0.0:
                events_per_nc_value = count / charge
                events_per_nC.append(events_per_nc_value)
                valid_runs.append(run_num)

                # Check if value is less than 0.08 and print run number
                if events_per_nc_value < 0.08:
                    print(f"Run {run_num}: {events_per_nc_value:.4f} events/nC")
                #endif
            else:
                print(f"Warning: Run {run_num} has zero charge, skipping.")
                zero_charge_count += 1
            #endif
        else:
            # Charge info missing from CSV for this run
            missing_charge_count += 1
        #endif
    #endfor

    # If no valid runs had charge info, skip plotting for this period
    if len(valid_runs) == 0:
        print(f"\nNo runs with charge information were found for period {period_label}.")
        if missing_charge_count > 0:
            print(f"  Note: {missing_charge_count} runs from this ROOT file are missing in the charge CSV.")
        #endif
        if zero_charge_count > 0:
            print(f"  Note: {zero_charge_count} runs had zero charge recorded.")
        #endif
        print("Skipping plot for this period.\n")
        #endif
        continue
    #endif

    # Sort runs and corresponding values by run number for a nicer plot
    valid_runs = np.array(valid_runs)
    events_per_nC = np.array(events_per_nC)
    sort_idx = np.argsort(valid_runs)
    valid_runs_sorted = valid_runs[sort_idx]
    events_per_nC_sorted = events_per_nC[sort_idx]

    # Create the plot with scattered points (no connecting line)
    plt.figure(figsize=(12, 6))
    plt.plot(valid_runs_sorted, events_per_nC_sorted, 'bo', markersize=4)
    plt.xlabel('Run Number')
    plt.ylabel('Events / nC')
    plt.title(f'DVCS Events per nC by Run Number ({period_label})')
    plt.grid(True, alpha=0.3)

    # Improve x-axis formatting for run numbers
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save the plot with a period-specific filename
    out_path = f'output/dvcs_per_nC_{period_label}.png'
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nPlot saved to {out_path}")
    print(f"Processed {len(valid_runs_sorted)} runs with valid charge information for {period_label}")
    if missing_charge_count > 0:
        print(f"{missing_charge_count} runs were missing from the charge CSV for this period.")
    #endif
    if zero_charge_count > 0:
        print(f"{zero_charge_count} runs had zero charge recorded and were skipped.")
    #endif
#endfof