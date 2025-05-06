import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

def generate_phase_space_plots(channel, correction, parent_dir, output_dir):
    """
    Generate a 1x3 canvas of 2D histograms for proton theta vs phi phase space.
    
    Parameters:
    - channel (str): 'ep' or 'eppi+pi-'
    - correction (str): Correction type ('noCorrections', 'timothy', 'mariana', 'krishna')
    - parent_dir (str): Directory containing input ROOT files
    - output_dir (str): Directory to save output PDFs
    """
    # Configuration based on channel
    if channel == 'ep':
        file_prefix = 'resIncl'
        theta_branch = 'p_theta'
        phi_branch = 'p_phi'
    elif channel == 'eppi+pi-':
        file_prefix = 'nSidis'
        theta_branch = 'p1_theta'
        phi_branch = 'p1_phi'
    else:
        raise ValueError(f"Invalid channel: {channel}. Must be 'ep' or 'eppi+pi-'")

    run_periods = {
        'fa18_inb': 'RGA Fa18 Inb',
        'fa18_out': 'RGA Fa18 Out',
        'sp19_inb': 'RGA Sp19 Inb'
    }

    # Create figure
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    plt.subplots_adjust(wspace=0.3, left=0.05, right=0.95)

    for idx, (run_key, run_label) in enumerate(run_periods.items()):
        # Construct file path
        filename = f"{file_prefix}_{run_key}_{channel}_{correction}.root"
        filepath = os.path.join(parent_dir, filename)
        
        if not os.path.exists(filepath):
            print(f"Warning: File not found - {filepath}")
            continue

        # Read data
        try:
            with uproot.open(filepath) as f:
                tree = f['PhysicsEvents']
                theta = np.degrees(tree[theta_branch].array(library='np'))
                phi = np.degrees(tree[phi_branch].array(library='np'))
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            continue

        # Create histogram
        ax = axs[idx]
        h = ax.hist2d(phi, theta, bins=100, range=[[-180, 180], [0, 180]], 
                     cmap='viridis', norm=LogNorm(), cmin=1)

        ax.set_xlabel(r'$\phi$ [deg]', fontsize=12, labelpad=10)
        ax.set_ylabel(r'$\theta$ [deg]', fontsize=12, labelpad=10)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_title(run_label, fontsize=14, pad=15)

        # Add colorbar to each subplot
        cb = fig.colorbar(h[3], ax=ax, pad=0.01)
        cb.set_label('Counts', fontsize=10)
        cb.ax.tick_params(labelsize=8)

    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{channel}_{correction}_phase_space.pdf")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

# Example usage
if __name__ == "__main__":
    PARENT_DIR = "/volatile/clas12/thayward/corrections_study/results/proton_energy_loss/"
    OUTPUT_DIR = "output/correction_study"

    # Analyze ep channel with different corrections
    for correction in ['noCorrections', 'timothy', 'mariana', 'krishna']:
        generate_phase_space_plots('ep', correction, PARENT_DIR, OUTPUT_DIR)

    # Analyze eppi+pi- channel with different corrections
    for correction in ['noCorrections', 'timothy', 'mariana', 'krishna']:
        generate_phase_space_plots('eppi+pi-', correction, PARENT_DIR, OUTPUT_DIR)