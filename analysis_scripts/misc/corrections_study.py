import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

def generate_phase_space_plots(channel, correction, plot_type, parent_dir, output_dir):
    """
    Generate a 1x3 canvas of 2D histograms for proton phase space analysis.
    
    Parameters:
    - channel (str): 'ep' or 'eppi+pi-'
    - correction (str): Correction type ('noCorrections', 'timothy', 'mariana', 'krishna')
    - plot_type (str): 'theta_phi' or 'theta_p'
    - parent_dir (str): Directory containing input ROOT files
    - output_dir (str): Directory to save output PDFs
    """
    # Configuration based on channel and plot type
    if channel == 'ep':
        file_prefix = 'resIncl'
        branch_config = {
            'theta_phi': ('p_theta', 'p_phi'),
            'theta_p': ('p_theta', 'p_p')
        }
    elif channel == 'eppi+pi-':
        file_prefix = 'nSidis'
        branch_config = {
            'theta_phi': ('p1_theta', 'p1_phi'),
            'theta_p': ('p1_theta', 'p1_p')
        }
    else:
        raise ValueError(f"Invalid channel: {channel}")

    try:
        y_branch, x_branch = branch_config[plot_type]
    except KeyError:
        raise ValueError(f"Invalid plot_type: {plot_type}. Must be 'theta_phi' or 'theta_p'")

    # Plot parameters
    plot_config = {
        'theta_phi': {
            'x_range': [0, 360],
            'y_range': [0, 140],
            'x_label': r'$\phi$ [deg]',
            'y_label': r'$\theta$ [deg]'
        },
        'theta_p': {
            'x_range': [0, 6],
            'y_range': [0, 140],
            'x_label': 'p [GeV]',
            'y_label': r'$\theta$ [deg]'
        }
    }

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

        # Read and process data
        try:
            with uproot.open(filepath) as f:
                tree = f['PhysicsEvents']
                
                # Get theta values and convert to degrees
                theta = np.degrees(tree[y_branch].array(library='np'))
                
                # Get x-axis values (CORRECTED LINE)
                x_data = tree[x_branch].array(library='np')  # Removed extra parenthesis
                
                # Special processing for phi
                if plot_type == 'theta_phi':
                    x_data = np.degrees(x_data) % 360  # Convert to degrees and wrap to 0-360

        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            continue

        # Create histogram
        ax = axs[idx]
        h = ax.hist2d(x_data, theta, bins=100,
                     range=[plot_config[plot_type]['x_range'], 
                            plot_config[plot_type]['y_range']],
                     cmap='viridis', norm=LogNorm(), cmin=1)

        ax.set_xlabel(plot_config[plot_type]['x_label'], fontsize=12, labelpad=10)
        ax.set_ylabel(plot_config[plot_type]['y_label'], fontsize=12, labelpad=10)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_title(run_label, fontsize=14, pad=15)

        # Add colorbar
        cb = fig.colorbar(h[3], ax=ax, pad=0.01)
        cb.set_label('Counts', fontsize=10)
        cb.ax.tick_params(labelsize=8)

    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 
                             f"{channel}_{correction}_{plot_type}_phase_space.pdf")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

# Main execution for noCorrections
if __name__ == "__main__":
    PARENT_DIR = "/volatile/clas12/thayward/corrections_study/results/proton_energy_loss/"
    OUTPUT_DIR = "output/correction_study"
    CORRECTION = 'noCorrections'

    # Generate both plot types for each channel
    for channel in ['ep', 'eppi+pi-']:
        for plot_type in ['theta_phi', 'theta_p']:
            generate_phase_space_plots(
                channel=channel,
                correction=CORRECTION,
                plot_type=plot_type,
                parent_dir=PARENT_DIR,
                output_dir=OUTPUT_DIR
            )