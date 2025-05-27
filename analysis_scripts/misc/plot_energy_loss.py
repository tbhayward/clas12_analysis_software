#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt

# === θ ranges ===
theta_fd = np.linspace(5, 39, 500)   # Forward Detector: 5–39°
theta_cd = np.linspace(25, 70, 500)  # Central Detector: 25–70°

# === Mariana ag matrix ===
mariana_theta_list = np.array([
    5.00, 7.00, 9.00, 11.00, 13.00,
    15.00, 17.00, 19.00, 21.00, 23.00,
    25.00, 27.00, 29.00, 31.00, 33.00
])
ag_matrix = np.array([
    [0.001039, -0.006952, -0.009509, -0.009879, -0.01279,
     -0.01157, -0.01018, -0.009222, -0.01355, -0.01207,
     -0.009474, -0.02216, -0.02105, -0.02118, -0.02360],
    [-0.0006922, 0.0009763, 0.001482, 0.001530, 0.002187,
     0.001953, 0.001688, 0.001668, 0.002849, 0.002495,
     0.001508, 0.004215, 0.003911, 0.003948, 0.004634],
    [0.0009806, 0.01157, 0.01485, 0.01588, 0.01945,
     0.01736, 0.01551, 0.01383, 0.01926, 0.01720,
     0.01464, 0.03250, 0.03231, 0.03296, 0.03608],
    [-0.008024, -0.01035, -0.01240, -0.01361, -0.01518,
     -0.01432, -0.01341, -0.01255, -0.01462, -0.01388,
     -0.01574, -0.02646, -0.02820, -0.03000, -0.03259]
])

def timothy_fd_fa18_inb(theta, p):
    A = 0.0099626 - 0.0002414*theta - 0.0000020*theta**2
    B = -0.01428267 + 0.00042833*theta + 0.00001081*theta**2
    C = 0.01197102 - 0.00055673*theta + 0.00000785*theta**2
    return A + B/p + C/p**2

def timothy_fd_fa18_out(theta, p):
    A = 0.0135790 - 0.0005303*theta
    B = -0.02165929 + 0.00121123*theta
    return A + B/p

def timothy_fd_sp19_inb(theta, p):
    A = 0.0095205 - 0.0001914*theta - 0.0000031*theta**2
    B = -0.01365658 + 0.00036322*theta + 0.00001217*theta**2
    C = 0.01175256 - 0.00053407*theta + 0.00000742*theta**2
    return A + B/p + C/p**2

def timothy_fd(theta, p):
    return timothy_fd_fa18_inb(theta, p)

def krishna_fd(theta, p):
    if theta < 27.0:
        if p < 2.4:
            return 0.001046*p**4 - 0.010446*p**3 + 0.036945*p**2 - 0.055368*p + 0.034539
        else:
            return 0.004741
    else:
        if p < 2.4:
            return 0.005519*p**4 - 0.046289*p**3 + 0.137504*p**2 - 0.177027*p + 0.094555
        else:
            return 0.004899

def mariana_fd(theta_vals, p):
    dp_table = -p * (
        ag_matrix[0] +
        ag_matrix[1]*p +
        ag_matrix[2]/p +
        ag_matrix[3]/(p**2)
    )
    idx = np.abs(theta_vals[:, None] - mariana_theta_list[None, :]).argmin(axis=1)
    return dp_table[idx]

def timothy_cd(theta, p):
    A = -0.2383991 + 0.0124992*theta - 0.0001646*theta**2
    B = 0.60123885 - 0.03128464*theta + 0.00041314*theta**2
    C = -0.44080146 + 0.02209857*theta - 0.00028224*theta**2
    return A + B*p + C*p**2

def stefan_pi_plus_fd(theta, p):
    if theta < 27 and p < 2.3:
        dp = 0.00389945 - 0.004062*p + 0.00321842*p**2 - 0.000698299*p**3 + 0*p**4
    elif theta < 27 and p >= 2.3:
        dp = 0.00389945 - 0.004062*2.3 + 0.00321842*2.3**2 - 0.000698299*2.3**3 + 0*2.3**4
    elif 27 < theta < 28 and p < 1.7:
        dp = 0.00727132 - 0.0117989*p + 0.00962999*p**2 - 0.00267005*p**3 + 0*p**4
    elif 27 < theta < 28 and p >= 1.7:
        dp = 0.00727132 - 0.0117989*1.7 + 0.00962999*1.7**2 - 0.00267005*1.7**3 + 0*1.7**4
    elif 28 < theta < 29 and p < 2:
        dp = 0.00844551 - 0.0128097*p + 0.00945956*p**2 - 0.00237992*p**3 + 0*p**4
    elif 28 < theta < 29 and p >= 2:
        dp = 0.00844551 - 0.0128097*2 + 0.00945956*2**2 - 0.00237992*2**3 + 0*2**4
    elif 29 < theta < 30 and p < 1.9:
        dp = 0.00959007 - 0.0139218*p + 0.0122966*p**2 - 0.0034012*p**3 + 0*p**4
    elif 29 < theta < 30 and p >= 1.9:
        dp = 0.00959007 - 0.0139218*1.9 + 0.0122966*1.9**2 - 0.0034012*1.9**3 + 0*1.9**4
    elif 30 < theta < 31 and p < 1.9:
        dp = 0.00542816 - 5.10739e-05*p + 0.000572038*p**2 - 0.000488883*p**3 + 0*p**4
    elif 30 < theta < 31 and p >= 1.9:
        dp = 0.00542816 - 5.10739e-05*1.9 + 0.000572038*1.9**2 - 0.000488883*1.9**3 + 0*1.9**4
    elif 31 < theta < 32 and p < 1.8:
        dp = 0.0060391 - 0.000516936*p - 0.00286595*p**2 + 0.00136604*p**3 + 0*p**4
    elif 31 < theta < 32 and p >= 1.8:
        dp = 0.0060391 - 0.000516936*1.8 - 0.00286595*1.8**2 + 0.00136604*1.8**3 + 0*1.8**4
    elif 32 < theta < 33 and p < 1.6:
        dp = 0.0140305 - 0.0285832*p + 0.0248799*p**2 - 0.00701311*p**3 + 0*1.6**4
    elif 32 < theta < 33 and p >= 1.6:
        dp = 0.0140305 - 0.0285832*1.6 + 0.0248799*1.6**2 - 0.00701311*1.6**3 + 0*1.6**4
    elif 33 < theta < 34 and p < 1.5:
        dp = 0.010815 - 0.0194244*p + 0.0174474*p**2 - 0.0049764*p**3 + 0*p**4
    elif 33 < theta < 34 and p >= 1.5:
        dp = 0.010815 - 0.0194244*1.5 + 0.0174474*1.5**2 - 0.0049764*1.5**3 + 0*1.5**4
    elif 34 < theta < 35 and p < 1.6:
        dp = 0.0105522 - 0.0176248*p + 0.0161142*p**2 - 0.00472288*p**3 + 0*p**4
    elif 34 < theta < 35 and p >= 1.6:
        dp = 0.0105522 - 0.0176248*1.6 + 0.0161142*1.6**2 - 0.00472288*1.6**3 + 0*1.6**4
    elif 35 < theta < 36 and p < 1.5:
        dp = 0.0103938 - 0.0164003*p + 0.0164045*p**2 - 0.00517012*p**3 + 0*p**4
    elif 35 < theta < 36 and p >= 1.5:
        dp = 0.0103938 - 0.0164003*1.5 + 0.0164045*1.5**2 - 0.00517012*1.5**3 + 0*1.5**4
    elif 36 < theta < 37 and p < 1.8:
        dp = (0.0441471 - 0.183937*p + 0.338784*p**2 - 0.298985*p**3 + 0.126905*p**4 - 0.0208286*p**5)
    elif 36 < theta < 37 and p >= 1.8:
        dp = (0.0441471 - 0.183937*1.8 + 0.338784*1.8**2 - 0.298985*1.8**3 + 0.126905*1.8**4 - 0.0208286*1.8**5)
    elif 37 < theta < 38 and p < 1.7:
        dp = (0.0726119 - 0.345004*p + 0.697789*p**2 - 0.685948*p**3 + 0.327195*p**4 - 0.0605621*p**5)
    elif 37 < theta < 38 and p >= 1.7:
        dp = (0.0726119 - 0.345004*1.7 + 0.697789*1.7**2 - 0.685948*1.7**3 + 0.327195*1.7**4 - 0.0605621*1.7**5)
    elif 38 < theta < 39 and p < 1.6:
        dp = 0.0247648 - 0.0797376*p + 0.126535*p**2 - 0.086545*p**3 + 0.0219304*p**4
    elif 38 < theta < 39 and p >= 1.6:
        dp = 0.0247648 - 0.0797376*1.6 + 0.126535*1.6**2 - 0.086545*1.6**3 + 0.0219304*1.6**4
    elif 39 < theta < 40 and p < 1.2:
        dp = 0.0208867 - 0.0492068*p + 0.0543187*p**2 - 0.0183393*p**3 + 0*p**4
    elif 39 < theta < 40 and p >= 1.2:
        dp = 0.0208867 - 0.0492068*1.2 + 0.0543187*1.2**2 - 0.0183393*1.2**3 + 0*1.2**4
    elif 40 < theta < 41 and p < 1.0:
        dp = 0.0148655 - 0.0203483*p + 0.00835867*p**2 + 0.00697134*p**3 + 0*p**4
    elif 40 < theta < 41 and p >= 1.0:
        dp = 0.0148655 - 0.0203483*1.0 + 0.00835867*1.0**2 + 0.00697134*1.0**3 + 0*1.0**4
    elif 41 < theta < 42 and p < 0.7:
        dp = 0.0223585 - 0.0365262*p - 0.0150027*p**2 + 0.0854164*p**3 - 0.0462718*p**4
    elif 41 < theta < 42 and p >= 0.7:
        dp = 0.007617
    elif theta > 42 and p < 0.75:
        dp = 0.0152373 - 0.0106377*p - 0.0257573*p**2 + 0.0344851*p**3 + 0*p**4
    elif theta > 42 and p >= 0.75:
        dp = 0.0152373 - 0.0106377*0.75 - 0.0257573*0.75**2 + 0.0344851*0.75**3 + 0*0.75**4
    else:
        dp = 0.0
    #endif
    return dp
#enddef

def krishna_pi_minus_fd(theta, p):
    if theta < 27.0:
        return p*0.00046571 + 0.00322164
    else:
        if p < 1.7:
            return -0.0024313*p**3 + 0.0094416*p**2 - 0.01257967*p + 0.0122432
        else:
            return 0.006199071

def main():
    out_dir = 'output'
    os.makedirs(out_dir, exist_ok=True)

    # Forward Detector proton plots
    fig_fd, axs_fd = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs_fd, [0.75, 1.75, 2.75]):
        ax.plot(theta_fd, timothy_fd(theta_fd, p), label='Timothy', linewidth=2)
        ax.plot(theta_fd, [krishna_fd(t, p) for t in theta_fd],
                label='Krishna', linestyle='--', linewidth=2)
        ax.plot(theta_fd, mariana_fd(theta_fd, p),
                label='Mariana', linestyle=':', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(5, 39); ax.set_ylim(-0.02, 0.03)
        ax.set_title(f'p = {p:.2f} GeV'); ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs_fd[0]: ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='lower left', frameon=True)
    fig_fd.suptitle('Forward Detector Proton Corrections')
    fig_fd.tight_layout(rect=[0,0,1,0.95]); fig_fd.savefig(f'{out_dir}/forward_detector.png')
    plt.close(fig_fd)

    # Central Detector proton plots
    fig_cd, axs_cd = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs_cd, [0.4, 0.75, 1.1]):
        ax.plot(theta_cd, timothy_cd(theta_cd, p), label='Timothy', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(25, 70); ax.set_ylim(-0.03, 0.03)
        ax.set_title(f'p = {p:.2f} GeV'); ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs_cd[0]: ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='upper right', frameon=True)
    fig_cd.suptitle('Central Detector Proton Corrections')
    fig_cd.tight_layout(rect=[0,0,1,0.95]); fig_cd.savefig(f'{out_dir}/central_detector.png')
    plt.close(fig_cd)

    # Timothy run-periods proton plots
    fig_rp, axs_rp = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs_rp, [0.75, 1.75, 2.75]):
        ax.plot(theta_fd, timothy_fd_fa18_inb(theta_fd, p), label='Fa18 Inb', linewidth=2)
        ax.plot(theta_fd, timothy_fd_fa18_out(theta_fd, p), label='Fa18 Out', linestyle='--', linewidth=2)
        ax.plot(theta_fd, timothy_fd_sp19_inb(theta_fd, p), label='Sp19 Inb', linestyle=':', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(5, 39); ax.set_ylim(-0.02, 0.02)
        ax.set_title(f'p = {p:.2f} GeV'); ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs_rp[0]: ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='lower left', frameon=True)
    fig_rp.suptitle("Timothy's FD Corrections Across Run Periods")
    fig_rp.tight_layout(rect=[0,0,1,0.95]); fig_rp.savefig(f'{out_dir}/timothy_run_periods.png')
    plt.close(fig_rp)

    # Pion comparison FD plots (p = 1,2,3 GeV)
    fig_pi, axs_pi = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs_pi, [1.0, 2.0, 3.0]):
        ax.plot(theta_fd, [stefan_pi_plus_fd(t, p) for t in theta_fd], label='Stefan π⁺ (Outb)', linewidth=2)
        ax.plot(theta_fd, [krishna_pi_minus_fd(t, p) for t in theta_fd], label='Krishna π⁻ (Inb)', linestyle='--', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(5, 39); ax.set_ylim(-0.01, 0.01)
        ax.set_title(f'p = {p:.1f} GeV'); ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs_pi[0]: ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='lower left', frameon=True)
    fig_pi.suptitle('π⁺ vs π⁻ Energy Loss (FD)')
    fig_pi.tight_layout(rect=[0,0,1,0.95]); fig_pi.savefig(f'{out_dir}/pi_comparison_fd.png')
    plt.close(fig_pi)

    print(f'Saved figures to {out_dir}/')

if __name__ == '__main__':
    main()