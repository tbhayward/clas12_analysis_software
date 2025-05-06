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


def plot_mx2_comparison(parent_dir, output_dir):
    """
    Analyzes missing mass squared (Mx²) spectra for different corrections
    """
    # Hardcoded configurations
    detectors = {
        1: {
            'name': 'Forward',
            'theta_bins': [0,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,80],
            'theta_labels': ['All θ'] + [
                '5-8', '8-11', '11-14', '14-17', '17-20', '20-23', '23-26',
                '26-29', '29-32', '32-35', '35-38', '38-41', '41-44', '44-47',
                '47-50', '50-53', '53-80'
            ]
        },
        2: {
            'name': 'Central',
            'theta_bins': [0,24,27,30,33,36,39,42,45,48,51,54,57,60,63,180],
            'theta_labels': ['All θ'] + [
                '24-27', '27-30', '30-33', '33-36', '36-39', '39-42', '42-45',
                '45-48', '48-51', '51-54', '54-57', '57-60', '60-63', '63-180'
            ]
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
    max_debug_events = 5
    
    # Mx² histogram parameters
    mx2_bins = np.linspace(0.3, 1.1, 100)
    
    for run in run_periods:
        for det_num, det_config in detectors.items():
            # Create figure with proper subplot grid
            n_total_plots = len(det_config['theta_labels'])
            n_cols = 4
            n_rows = int(np.ceil((n_total_plots-1)/n_cols)) + 1
            
            fig = plt.figure(figsize=(20, 5*n_rows))
            gs = gridspec.GridSpec(n_rows, n_cols, figure=fig)
            
            # Load all data first
            all_data = {}
            for corr in corrections:
                filename = f"resIncl_{run}_ep_{corr}.root"
                filepath = os.path.join(parent_dir, filename)
                
                try:
                    with uproot.open(filepath) as f:
                        tree = f['PhysicsEvents']
                        data = tree.arrays(['Mx2', 'p_p', 'p_theta', 'detector'], library='np')
                        mask = (data['detector'] == det_num)
                        
                        if np.sum(mask) == 0:
                            print(f"No events in {filename} for {det_config['name']}")
                            continue
                            
                        all_data[corr] = {
                            'Mx2': data['Mx2'][mask],
                            'p_p': data['p_p'][mask],
                            'theta': np.degrees(data['p_theta'][mask])
                        }
                        
                except Exception as e:
                    print(f"Error loading {filepath}: {str(e)}")
                    continue

            # Debug print proton momenta AND Mx² values
            if all_data:
                print(f"\n{' Event Data Comparison ':=^80}")
                print(f"{'Event':5} | {'p_p (GeV)':^50} | {'Mx² (GeV²)':^50}")
                print("-"*130)
                
                valid_corrections = [corr for corr in corrections if corr in all_data]
                min_events = min([len(all_data[corr]['p_p']) for corr in valid_corrections]) if valid_corrections else 0
                n_print = min(max_debug_events, min_events)
                
                for i in range(n_print):
                    # Print p_p
                    line_p = f"{i:5} | "
                    for corr in corrections:
                        if corr in all_data and i < len(all_data[corr]['p_p']):
                            line_p += f"{corr}: {all_data[corr]['p_p'][i]:.3f} | "
                        else:
                            line_p += f"{corr}: N/A | "
                    
                    # Print Mx²
                    line_mx2 = f"{' ':<5} | "
                    for corr in corrections:
                        if corr in all_data and i < len(all_data[corr]['Mx2']):
                            line_mx2 += f"{corr}: {all_data[corr]['Mx2'][i]:.3f} | "
                        else:
                            line_mx2 += f"{corr}: N/A | "
                    
                    print(line_p)
                    print(line_mx2)
                    print("-"*130)

            # Plot integrated spectrum
            ax_int = fig.add_subplot(gs[0, :])
            for (corr, color, ls) in zip(corrections, colors, line_styles):
                if corr not in all_data:
                    continue
                ax_int.hist(
                    all_data[corr]['Mx2'], bins=mx2_bins,
                    histtype='step', color=color, linestyle=ls, linewidth=2,
                    label=corr_labels[corr], density=False
                )
            
            ax_int.set(xlabel=r'$M_{x}^{2}$ (GeV²)', ylabel='Counts',
                      title=f"{det_config['name']} Detector - {run} - Integrated",
                      xlim=(0.3, 1.1))
            ax_int.legend()
            ax_int.grid(True, alpha=0.3)
            
            # Plot theta-binned spectra
            for idx in range(1, n_total_plots):
                row = (idx-1) // n_cols + 1
                col = (idx-1) % n_cols
                ax = fig.add_subplot(gs[row, col])
                
                theta_min = det_config['theta_bins'][idx]
                theta_max = det_config['theta_bins'][idx+1]
                
                # Plot all corrections
                artists = []
                for (corr, color, ls) in zip(corrections, colors, line_styles):
                    if corr not in all_data:
                        continue
                    
                    mask = (all_data[corr]['theta'] >= theta_min) & (all_data[corr]['theta'] < theta_max)
                    mx2_data = all_data[corr]['Mx2'][mask]
                    
                    if len(mx2_data) > 0:
                        n, bins, patches = ax.hist(
                            mx2_data, bins=mx2_bins,
                            histtype='step', color=color, linestyle=ls, linewidth=2,
                            label=corr_labels[corr], density=False
                        )
                        artists.append(patches[0])
                
                ax.set(xlabel=r'$M_{x}^{2}$ (GeV²)', ylabel='Counts',
                      title=f'θ: {det_config["theta_labels"][idx]}°',
                      xlim=(0.3, 1.1))
                ax.grid(True, alpha=0.3)
                
                if artists:
                    ax.legend(handles=artists, fontsize=8)

            plt.tight_layout()
            output_file = os.path.join(output_dir, f'Mx2_comparison_{run}_{det_config["name"]}.pdf')
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

    # Generate Mx2 comparison plots
    plot_mx2_comparison(PARENT_DIR, OUTPUT_DIR)