import uproot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Create output directory if it doesn't exist
os.makedirs('output', exist_ok=True)

# Read the CSV file with run information
csv_file = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Read CSV, skipping comment lines
run_info_df = pd.read_csv(csv_file, comment='#', header=None, 
                         names=['runnum', 'charge_nC', 'col2', 'col3', 'col4', 'col5'])

# Create a dictionary mapping run number to accumulated charge
run_charge_map = dict(zip(run_info_df['runnum'], run_info_df['charge_nC']))

# Open the ROOT file and access the PhysicsEvents tree
root_file = uproot.open("/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root")
tree = root_file["PhysicsEvents"]

# Get the runnum branch data
runnum_branch = tree["runnum"]
runnum_data = runnum_branch.array()

# Count events per run number
unique_runs, event_counts = np.unique(runnum_data, return_counts=True)

# Calculate events per nC for each run
events_per_nC = []
valid_runs = []

print("Runs with events/nC < 0.8:")
print("-" * 30)

for run_num, count in zip(unique_runs, event_counts):
    if run_num in run_charge_map:
        charge = run_charge_map[run_num]
        if charge > 0:  # Avoid division by zero
            events_per_nc_value = count / charge
            events_per_nC.append(events_per_nc_value)
            valid_runs.append(run_num)
            
            # Check if value is less than 0.8 and print run number
            if events_per_nc_value < 0.8:
                print(f"Run {run_num}: {events_per_nc_value:.4f} events/nC")
        else:
            print(f"Warning: Run {run_num} has zero charge, skipping")
    else:
        print(f"Warning: Run {run_num} not found in CSV file, skipping")

# Create the plot with scattered points (no connecting line)
plt.figure(figsize=(12, 6))
plt.plot(valid_runs, events_per_nC, 'bo', markersize=4)  # Removed the '-'' to remove connecting line
plt.xlabel('Run Number')
plt.ylabel('Events / nC')
plt.title('DVCS Events per nC by Run Number')
plt.grid(True, alpha=0.3)

# Improve x-axis formatting for run numbers
plt.xticks(rotation=45)
plt.tight_layout()

# Save the plot
plt.savefig('output/dvcs_per_nC.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"\nPlot saved to output/dvcs_per_nC.png")
print(f"Processed {len(valid_runs)} runs with valid charge information")