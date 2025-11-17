import uproot
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import colors

# Create output directory
os.makedirs('output/DC_HV_scan', exist_ok=True)

# List of root files - fixed the typo in the path
file_paths = [
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_10_10.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_10_11.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_11_10.root',  # Fixed path
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_11_11.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_11_12.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_12_10.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_12_11.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_11_11_11.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_12_13_13.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_9_10_10.root'
]

# Initialize arrays to store all data
all_Mx = []
all_Q2 = []
all_x = []
all_e_theta = []
all_e_p = []
all_p_theta = []
all_p_p = []

# Loop over all files and read trees
for file_path in file_paths:
    try:
        with uproot.open(file_path) as file:
            tree = file['PhysicsEvents']
            
            # Read the required branches - using Mx2 instead of Mx
            Mx2 = tree['Mx2'].array()
            Q2 = tree['Q2'].array()
            x = tree['x'].array()
            e_theta = tree['e_theta'].array()
            e_p = tree['e_p'].array()
            p_theta = tree['p_theta'].array()
            p_p = tree['p_p'].array()
            
            # Calculate Mx from Mx2 (taking square root)
            Mx = np.sqrt(Mx2)
            
            # Apply Mx cut (0.7 to 1.1 GeV)
            mask = (Mx >= 0.7) & (Mx <= 1.1)
            
            # Append data after applying cut
            all_Mx.extend(Mx[mask])
            all_Q2.extend(Q2[mask])
            all_x.extend(x[mask])
            all_e_theta.extend(e_theta[mask])
            all_e_p.extend(e_p[mask])
            all_p_theta.extend(p_theta[mask])
            all_p_p.extend(p_p[mask])
            
            print(f"Processed {file_path}: {np.sum(mask)} events after Mx cut")
            
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

# Convert to numpy arrays for easier handling
all_Mx = np.array(all_Mx)
all_Q2 = np.array(all_Q2)
all_x = np.array(all_x)
all_e_theta = np.array(all_e_theta)
all_e_p = np.array(all_e_p)
all_p_theta = np.array(all_p_theta)
all_p_p = np.array(all_p_p)

print(f"Total events after Mx cut: {len(all_Mx)}")

# Check if we have data before plotting
if len(all_Mx) == 0:
    print("No events found after cuts! Check your data and cuts.")
    exit()

# Create the 2x2 plot
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Top left: Mx histogram
axes[0,0].hist(all_Mx, bins=100, range=(0.7, 1.1), color='blue', alpha=0.7, edgecolor='black')
axes[0,0].set_xlabel('M$_x$ (GeV)')
axes[0,0].set_ylabel('Counts')
axes[0,0].grid(True, alpha=0.3)
axes[0,0].set_title('Missing Mass Distribution')

# Top right: Q2 vs x
if len(all_x) > 0 and len(all_Q2) > 0:
    hist1 = axes[0,1].hist2d(all_x, all_Q2, bins=100, range=[[0, 0.8], [1, 9]], 
                             cmap='viridis', norm=colors.LogNorm())
    axes[0,1].set_xlabel('x$_B$')
    axes[0,1].set_ylabel('Q$^{2}$ (GeV$^{2}$)')
    axes[0,1].grid(True, alpha=0.3)
    axes[0,1].set_title('Q$^{2}$ vs x$_B$')
    plt.colorbar(hist1[3], ax=axes[0,1])

# Bottom left: e_theta vs e_p
if len(all_e_p) > 0 and len(all_e_theta) > 0:
    hist2 = axes[1,0].hist2d(all_e_p, all_e_theta, bins=100, cmap='viridis', norm=colors.LogNorm())
    axes[1,0].set_xlabel('e$_p$ (GeV)')
    axes[1,0].set_ylabel('e$_{\\theta}$')
    axes[1,0].grid(True, alpha=0.3)
    axes[1,0].set_title('Electron $\\theta$ vs Electron Momentum')
    plt.colorbar(hist2[3], ax=axes[1,0])

# Bottom right: p_theta vs p_p - fixed escape sequences
if len(all_p_p) > 0 and len(all_p_theta) > 0:
    hist3 = axes[1,1].hist2d(all_p_p, all_p_theta, bins=100, cmap='viridis', norm=colors.LogNorm())
    axes[1,1].set_xlabel(r'$\pi^{+}_{p}$ (GeV)')  # raw string for LaTeX
    axes[1,1].set_ylabel(r'$\pi^{+}_{\theta}$')   # raw string for LaTeX
    axes[1,1].grid(True, alpha=0.3)
    axes[1,1].set_title(r'$\pi^{+}$ $\theta$ vs $\pi^{+}$ Momentum')  # raw string for LaTeX
    plt.colorbar(hist3[3], ax=axes[1,1])

# Adjust layout and save
plt.tight_layout()
plt.savefig('output/DC_HV_scan/kinematics.png', dpi=300, bbox_inches='tight')
plt.close()

print("Plots saved to output/DC_HV_scan/kinematics.png")