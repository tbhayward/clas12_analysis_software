#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# === θ ranges ===
theta_fd = np.linspace(5, 33, 500)   # Forward Detector: 5–33°
theta_cd = np.linspace(25, 70, 500)  # Central Detector: 25–70°

# === Mariana ag matrix ===
mariana_theta_list = np.array([5.00, 7.00, 9.00, 11.00, 13.00, 15.00,
                               17.00, 19.00, 21.00, 23.00, 25.00,
                               27.00, 29.00, 31.00, 33.00])
ag_matrix = np.array([
    [0.001039, -0.006952, -0.009509, -0.009879, -0.01279, -0.01157, -0.01018,
     -0.009222, -0.01355, -0.01207, -0.009474, -0.02216, -0.02105,
     -0.02118, -0.02360],
    [-0.0006922, 0.0009763, 0.001482, 0.001530, 0.002187, 0.001953,
     0.001688, 0.001668, 0.002849, 0.002495, 0.001508, 0.004215,
     0.003911, 0.003948, 0.004634],
    [0.0009806, 0.01157, 0.01485, 0.01588, 0.01945, 0.01736,
     0.01551, 0.01383, 0.01926, 0.01720, 0.01464, 0.03250,
     0.03231, 0.03296, 0.03608],
    [-0.008024, -0.01035, -0.01240, -0.01361, -0.01518, -0.01432,
     -0.01341, -0.01255, -0.01462, -0.01388, -0.01574, -0.02646,
     -0.02820, -0.03000, -0.03259]
])

# === Correction functions ===
def timothy_fd(theta, p):
    A = 0.0099626 - 0.0002414*theta - 0.0000020*theta**2
    B = -0.01428267 + 0.00042833*theta + 0.00001081*theta**2
    C = 0.01197102 - 0.00055673*theta + 0.00000785*theta**2
    return A + B/p + C/p**2

def krishna_fd(theta, p):
    if theta < 27:
        return 0.001046*p**4 - 0.010446*p**3 + 0.036945*p**2 - 0.055368*p + 0.034539
    else:
        return 0.005519*p**4 - 0.046289*p**3 + 0.137504*p**2 - 0.177027*p + 0.094555

def mariana_fd(theta_vals, p):
    dp_vals = -p*(ag_matrix[0] + ag_matrix[1]*p + ag_matrix[2]/p + ag_matrix[3]/(p**2))
    f = interp1d(mariana_theta_list, dp_vals,
                 kind='quadratic', bounds_error=False,
                 fill_value='extrapolate')
    return f(theta_vals)

def timothy_cd(theta, p):
    A = -0.2383991 + 0.0124992*theta - 0.0001646*theta**2
    B = 0.60123885 - 0.03128464*theta + 0.00041314*theta**2
    C = -0.44080146 + 0.02209857*theta - 0.00028224*theta**2
    return A + B*p + C*p**2

def main():
    out_dir = 'output'
    os.makedirs(out_dir, exist_ok=True)

    # === Forward Detector ===
    fig_fd, axs_fd = plt.subplots(1, 3, figsize=(12, 4),
                                  sharey=True,
                                  gridspec_kw={'wspace': 0})
    for ax, p in zip(axs_fd, [1.0, 2.0, 3.0]):
        ax.plot(theta_fd, timothy_fd(theta_fd, p),
                label='Timothy', linewidth=2)
        ax.plot(theta_fd, [krishna_fd(t, p) for t in theta_fd],
                label='Krishna', linewidth=2, linestyle='--')
        ax.plot(theta_fd, mariana_fd(theta_fd, p),
                label='Mariana', linewidth=2, linestyle=':')
        ax.set_xlim(5, 33)
        ax.set_ylim(-0.01, 0.03)
        ax.set_title(f'p = {p:.1f} GeV', fontsize=12)
        ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs_fd[0]:
            ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='bottom right', frameon=True)

    fig_fd.suptitle('Forward Detector Proton Energy Loss Corrections', fontsize=14)
    fig_fd.tight_layout(rect=[0, 0, 1, 0.95])
    fig_fd.savefig(f'{out_dir}/forward_detector.png')

    # === Central Detector ===
    fig_cd, axs_cd = plt.subplots(1, 3, figsize=(12, 4),
                                  sharey=True,
                                  gridspec_kw={'wspace': 0})
    for ax, p in zip(axs_cd, [0.4, 0.75, 1.1]):
        ax.plot(theta_cd, timothy_cd(theta_cd, p),
                label='Timothy', linewidth=2)
        ax.set_xlim(25, 70)
        ax.set_ylim(-0.02, 0.02)
        ax.set_title(f'p = {p:.2f} GeV', fontsize=12)
        ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs_cd[0]:
            ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='upper right', frameon=True)

    fig_cd.suptitle('Central Detector Proton Energy Loss Corrections', fontsize=14)
    fig_cd.tight_layout(rect=[0, 0, 1, 0.95])
    fig_cd.savefig(f'{out_dir}/central_detector.png')

    print(f'Saved figures to {out_dir}/')

if __name__ == '__main__':
    main()