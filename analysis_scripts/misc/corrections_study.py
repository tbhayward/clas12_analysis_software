import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import os

def generate_phase_space_plots(channel, correction, plot_type, parent_dir, output_dir):
    """
    Generate a 1x3 canvas of 2D histograms for proton phase space analysis.
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
        filename = f"{file_prefix}_{run_key}_{channel}_{correction}.root"
        filepath = os.path.join(parent_dir, filename)
        
        if not os.path.exists(filepath):
            print(f"Warning: File not found - {filepath}")
            continue

        try:
            with uproot.open(filepath) as f:
                tree = f['PhysicsEvents']
                theta = np.degrees(tree[y_branch].array(library='np'))
                x_data = tree[x_branch].array(library='np')
                
                if plot_type == 'theta_phi':
                    x_data = np.degrees(x_data) % 360

        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            continue

        ax = axs[idx]
        h = ax.hist2d(x_data, theta, bins=100,
                     range=[plot_config[plot_type]['x_range'], 
                            plot_config[plot_type]['y_range']],
                     cmap='viridis', norm=LogNorm(), cmin=1)

        ax.set_xlabel(plot_config[plot_type]['x_label'], fontsize=12)
        ax.set_ylabel(plot_config[plot_type]['y_label'], fontsize=12)
        ax.set_title(run_label, fontsize=14)
        fig.colorbar(h[3], ax=ax).set_label('Counts', fontsize=10)

    output_file = os.path.join(output_dir, f"{channel}_{correction}_{plot_type}_phase_space.pdf")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

def plot_mx2_comparison(parent_dir, output_dir):
    """
    Analyzes missing mass squared (Mx²) spectra for different corrections
    """
    detectors = {
        1: {
            'name': 'Forward',
            'theta_bins': [0, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 80]
        },
        2: {
            'name': 'Central',
            'theta_bins': [0, 36, 39, 42, 45, 48, 51, 54, 57, 180]
        }
    }

    # build theta_labels for each detector, skipping the first (0–edge) bin
    for det_config in detectors.values():
        bins = det_config['theta_bins']
        det_config['theta_labels'] = ['All θ'] + [
            f"{bins[i]}-{bins[i+1]}" for i in range(1, len(bins)-1)
        ]
    #endfor

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

    mx2_bins = np.linspace(-0.4, 0.4, 100)

    for run in run_periods:
        for det_num, det_config in detectors.items():
            n_total_plots = len(det_config['theta_labels'])
            n_cols = 4
            n_rows = int(np.ceil((n_total_plots - 1) / n_cols)) + 1

            fig = plt.figure(figsize=(20, 5 * n_rows))
            gs = gridspec.GridSpec(n_rows, n_cols, figure=fig)

            all_data = {}
            for corr in corrections:
                filename = f"resIncl_{run}_ep_{corr}.root"
                filepath = os.path.join(parent_dir, filename)

                try:
                    with uproot.open(filepath) as f:
                        tree = f['PhysicsEvents']
                        data = tree.arrays(['Mx2', 'p_theta', 'detector'], library='np')
                        mask = (data['detector'] == det_num)

                        if mask.sum() > 0:
                            all_data[corr] = {
                                'Mx2': data['Mx2'][mask],
                                'theta': np.degrees(data['p_theta'][mask])
                            }
                        #endif
                except Exception as e:
                    print(f"Error loading {filepath}: {e}")
                #endtry
            #endfor

            # Integrated spectrum plot
            ax_int = fig.add_subplot(gs[0, :])
            for corr, color, ls in zip(corrections, colors, line_styles):
                if corr in all_data:
                    ax_int.hist(
                        all_data[corr]['Mx2'], bins=mx2_bins,
                        histtype='step', color=color, linestyle=ls,
                        label=corr_labels[corr]
                    )
                #endif
            #endfor

            ax_int.set(
                xlabel=r'$M_{x}^{2}$ (GeV²)', ylabel='Counts',
                xlim=(-0.4, 0.4),
                title=f"{det_config['name']} Detector - {run}"
            )
            ax_int.legend()
            ax_int.grid(True, alpha=0.3)

            # Theta-binned spectra
            for idx in range(1, n_total_plots):
                row = (idx - 1) // n_cols + 1
                col = (idx - 1) % n_cols
                ax = fig.add_subplot(gs[row, col])

                theta_min = det_config['theta_bins'][idx]
                theta_max = det_config['theta_bins'][idx + 1]

                artists = []
                for corr, color, ls in zip(corrections, colors, line_styles):
                    if corr in all_data:
                        mask = (
                            (all_data[corr]['theta'] >= theta_min) &
                            (all_data[corr]['theta'] < theta_max)
                        )
                        mx2_data = all_data[corr]['Mx2'][mask]
                        if len(mx2_data) > 0:
                            _, _, patches = ax.hist(
                                mx2_data, bins=mx2_bins,
                                histtype='step', color=color, linestyle=ls,
                                label=corr_labels[corr]
                            )
                            artists.append(patches[0])
                        #endif
                    #endif
                #endfor

                ax.set(
                    xlabel=r'$M_{x}^{2}$ (GeV²)',
                    ylabel='Counts',
                    xlim=(-0.4, 0.4),
                    title=f'θ: {det_config["theta_labels"][idx]}°'
                )
                if artists:
                    ax.legend(handles=artists, fontsize=8)
                #endif
            #endfor

            output_file = os.path.join(
                output_dir,
                f'Mx2_comparison_{run}_{det_config["name"]}.pdf'
            )
            plt.savefig(output_file, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_file}")
        #endfor
    #endfor

if __name__ == "__main__":
    PARENT_DIR = "/volatile/clas12/thayward/corrections_study/results/proton_energy_loss/"
    OUTPUT_DIR = "output/correction_study"
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Generate phase space plots (uncomment to run)
    # for channel in ['ep', 'eppi+pi-']:
    #     for plot_type in ['theta_phi', 'theta_p']:
    #         generate_phase_space_plots(
    #             channel=channel,
    #             correction='noCorrections',
    #             plot_type=plot_type,
    #             parent_dir=PARENT_DIR,
    #             output_dir=OUTPUT_DIR
    #         )

    # Generate Mx² comparison plots
    try:
        # Generate plots
        plot_mx2_comparison(PARENT_DIR, OUTPUT_DIR)
        
        # Uncomment to run phase space plots
        # generate_phase_space_plots(...)
        
    finally:
        # Clean up matplotlib resources
        plt.close('all')