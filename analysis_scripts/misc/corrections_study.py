import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
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
            'x_range': [0, 7],
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


def plot_w_comparison(parent_dir, output_dir):
    """
    Plot W spectrum comparisons for different corrections and detectors
    """
    # Hardcoded bin configurations
    detectors = {
        1: {
            'name': 'Forward',
            'theta_bins': [0, 5] + list(range(5, 54, 3)) + [90],
            'theta_labels': ['0-5'] + [f'{x}-{x+3}' for x in range(5, 52, 3)] + ['53-90']
        },
        2: {
            'name': 'Central',
            'theta_bins': [0, 24] + list(range(24, 64, 3)) + [180],
            'theta_labels': ['0-24'] + [f'{x}-{x+3}' for x in range(24, 61, 3)] + ['63-180']
        }
    }
    
    corrections = ['noCorrections', 'timothy', 'krishna', 'mariana']
    corr_labels = {
        'noCorrections': 'No Corrections',
        'timothy': "Timothy's",
        'krishna': "Krishna's",
        'mariana': "Mariana's"
    }
    colors = ['black', 'red', 'green', 'blue']
    line_styles = ['-', '--', ':', '-.']
    run_periods = ['fa18_inb', 'fa18_out', 'sp19_inb']
    
    # W histogram parameters
    w_bins = np.linspace(0.6, 1.2, 100)
    
    for run in run_periods:
        for det_num, det_config in detectors.items():
            theta_bins = det_config['theta_bins']
            n_theta_bins = len(theta_bins) - 1
            
            # Create figure layout
            fig = plt.figure(figsize=(20, 16))
            gs = gridspec.GridSpec(5, 4, figure=fig)
            ax_main = fig.add_subplot(gs[0, 0:2])
            ax_theta = [fig.add_subplot(gs[i//2+1, i%2]) for i in range(8)]
            
            # Load all correction data
            all_data = {}
            for corr in corrections:
                filename = f"resIncl_{run}_ep_{corr}.root"
                filepath = os.path.join(parent_dir, filename)
                
                try:
                    with uproot.open(filepath) as f:
                        tree = f['PhysicsEvents']
                        data = tree.arrays(['W', 'p_theta', 'detector'], library='np')
                        mask = (data['detector'] == det_num)
                        
                        if np.sum(mask) == 0:
                            print(f"No events in {filename} for detector {det_num}")
                            continue
                            
                        all_data[corr] = {
                            'W': data['W'][mask],
                            'theta': np.degrees(data['p_theta'][mask])
                        }
                except Exception as e:
                    print(f"Error loading {filepath}: {str(e)}")
                    continue
            
            # Plot integrated spectrum with all corrections
            for corr, color, ls in zip(corrections, colors, line_styles):
                if corr not in all_data:
                    continue
                ax_main.hist(all_data[corr]['W'], bins=w_bins, histtype='step',
                            color=color, linestyle=ls, linewidth=1.5,
                            label=corr_labels[corr], density=False)
            
            ax_main.set(xlabel='W (GeV)', ylabel='Counts',
                      title=f'{det_config["name"]} Detector - {run}', xlim=(0.6, 1.2))
            ax_main.legend(prop={'size': 10})
            ax_main.grid(True, alpha=0.3)
            
            # Plot theta-binned spectra
            for idx in range(n_theta_bins):
                ax = ax_theta[idx]
                theta_min = theta_bins[idx]
                theta_max = theta_bins[idx+1]
                
                for corr, color, ls in zip(corrections, colors, line_styles):
                    if corr not in all_data:
                        continue
                        
                    mask = (all_data[corr]['theta'] >= theta_min) & (all_data[corr]['theta'] < theta_max)
                    w_data = all_data[corr]['W'][mask]
                    
                    if len(w_data) > 0:
                        ax.hist(w_data, bins=w_bins, histtype='step',
                               color=color, linestyle=ls, linewidth=1.5,
                               label=corr_labels[corr], density=False)
                
                ax.set(xlabel='W (GeV)', ylabel='Counts',
                      title=f'θ: {det_config["theta_labels"][idx]}°',
                      xlim=(0.6, 1.2))
                ax.grid(True, alpha=0.3)
                ax.tick_params(labelsize=9)
            
            plt.tight_layout()
            output_file = os.path.join(output_dir, f'W_comparison_{run}_{det_config["name"]}.pdf')
            plt.savefig(output_file, bbox_inches='tight')
            plt.close()
            print(f"Saved: {output_file}")

# Main execution for noCorrections
if __name__ == "__main__":
    PARENT_DIR = "/volatile/clas12/thayward/corrections_study/results/proton_energy_loss/"
    OUTPUT_DIR = "output/correction_study"
    CORRECTION = 'noCorrections'

    # # Generate both plot types for each channel
    # for channel in ['ep', 'eppi+pi-']:
    #     for plot_type in ['theta_phi', 'theta_p']:
    #         generate_phase_space_plots(
    #             channel=channel,
    #             correction=CORRECTION,
    #             plot_type=plot_type,
    #             parent_dir=PARENT_DIR,
    #             output_dir=OUTPUT_DIR
    #         )

    # Generate W comparison plots
    plot_w_comparison(PARENT_DIR, OUTPUT_DIR)