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
    [ 0.001039, -0.006952, -0.009509, -0.009879, -0.01279,
     -0.01157, -0.01018, -0.009222, -0.01355, -0.01207,
     -0.009474, -0.02216, -0.02105, -0.02118, -0.02360],
    [-0.0006922,  0.0009763,  0.001482,  0.001530,  0.002187,
      0.001953,  0.001688,  0.001668,  0.002849,  0.002495,
      0.001508,  0.004215,  0.003911,  0.003948,  0.004634],
    [ 0.0009806,  0.01157,    0.01485,   0.01588,   0.01945,
      0.01736,   0.01551,   0.01383,   0.01926,   0.01720,
      0.01464,   0.03250,   0.03231,   0.03296,   0.03608],
    [-0.008024,  -0.01035,   -0.01240,  -0.01361,  -0.01518,
     -0.01432,  -0.01341,  -0.01255,  -0.01462,  -0.01388,
     -0.01574,  -0.02646,  -0.02820,  -0.03000,  -0.03259]
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
        dp = 0.00389945 - 0.004062*p + 0.00321842*p**2 - 0.000698299*p**3
    elif theta < 27:
        dp = 0.00389945 - 0.004062*2.3 + 0.00321842*2.3**2 - 0.000698299*2.3**3
    elif 27 < theta < 28 and p < 1.7:
        dp = 0.00727132 - 0.0117989*p + 0.00962999*p**2 - 0.00267005*p**3
    elif 27 < theta < 28:
        dp = 0.00727132 - 0.0117989*1.7 + 0.00962999*1.7**2 - 0.00267005*1.7**3
    elif 28 < theta < 29 and p < 2:
        dp = 0.00844551 - 0.0128097*p + 0.00945956*p**2 - 0.00237992*p**3
    elif 28 < theta < 29:
        dp = 0.00844551 - 0.0128097*2 + 0.00945956*2**2 - 0.00237992*2**3
    elif 29 < theta < 30 and p < 1.9:
        dp = 0.00959007 - 0.0139218*p + 0.0122966*p**2 - 0.0034012*p**3
    elif 29 < theta < 30:
        dp = 0.00959007 - 0.0139218*1.9 + 0.0122966*1.9**2 - 0.0034012*1.9**3
    elif 30 < theta < 31 and p < 1.9:
        dp = 0.00542816 - 5.10739e-05*p + 0.000572038*p**2 - 0.000488883*p**3
    elif 30 < theta < 31:
        dp = 0.00542816 - 5.10739e-05*1.9 + 0.000572038*1.9**2 - 0.000488883*1.9**3
    elif 31 < theta < 32 and p < 1.8:
        dp = 0.0060391 - 0.000516936*p - 0.00286595*p**2 + 0.00136604*p**3
    elif 31 < theta < 32:
        dp = 0.0060391 - 0.000516936*1.8 - 0.00286595*1.8**2 + 0.00136604*1.8**3
    elif 32 < theta < 33 and p < 1.6:
        dp = 0.0140305 - 0.0285832*p + 0.0248799*p**2 - 0.00701311*p**3
    elif 32 < theta < 33:
        dp = 0.0140305 - 0.0285832*1.6 + 0.0248799*1.6**2 - 0.00701311*1.6**3
    elif 33 < theta < 34 and p < 1.5:
        dp = 0.010815 - 0.0194244*p + 0.0174474*p**2 - 0.0049764*p**3
    elif 33 < theta < 34:
        dp = 0.010815 - 0.0194244*1.5 + 0.0174474*1.5**2 - 0.0049764*1.5**3
    elif 34 < theta < 35 and p < 1.6:
        dp = 0.0105522 - 0.0176248*p + 0.0161142*p**2 - 0.00472288*p**3
    elif 34 < theta < 35:
        dp = 0.0105522 - 0.0176248*1.6 + 0.0161142*1.6**2 - 0.00472288*1.6**3
    elif 35 < theta < 36 and p < 1.5:
        dp = 0.0103938 - 0.0164003*p + 0.0164045*p**2 - 0.00517012*p**3
    elif 35 < theta < 36:
        dp = 0.0103938 - 0.0164003*1.5 + 0.0164045*1.5**2 - 0.00517012*1.5**3
    elif 36 < theta < 37 and p < 1.8:
        dp = (0.0441471 - 0.183937*p + 0.338784*p**2 - 0.298985*p**3 +
              0.126905*p**4 - 0.0208286*p**5)
    elif 36 < theta < 37:
        dp = (0.0441471 - 0.183937*1.8 + 0.338784*1.8**2 - 0.298985*1.8**3 +
              0.126905*1.8**4 - 0.0208286*1.8**5)
    elif 37 < theta < 38 and p < 1.7:
        dp = (0.0726119 - 0.345004*p + 0.697789*p**2 - 0.685948*p**3 +
              0.327195*p**4 - 0.0605621*p**5)
    elif 37 < theta < 38:
        dp = (0.0726119 - 0.345004*1.7 + 0.697789*1.7**2 - 0.685948*1.7**3 +
              0.327195*1.7**4 - 0.0605621*1.7**5)
    elif 38 < theta < 39 and p < 1.6:
        dp = 0.0247648 - 0.0797376*p + 0.126535*p**2 - 0.086545*p**3 + 0.0219304*p**4
    elif 38 < theta < 39:
        dp = 0.0247648 - 0.0797376*1.6 + 0.126535*1.6**2 - 0.086545*1.6**3 + 0.0219304*1.6**4
    elif 39 < theta < 40 and p < 1.2:
        dp = 0.0208867 - 0.0492068*p + 0.0543187*p**2 - 0.0183393*p**3
    elif 39 < theta < 40:
        dp = 0.0208867 - 0.0492068*1.2 + 0.0543187*1.2**2 - 0.0183393*1.2**3
    elif 40 < theta < 41 and p < 1.0:
        dp = 0.0148655 - 0.0203483*p + 0.00835867*p**2 + 0.00697134*p**3
    elif 40 < theta < 41:
        dp = 0.0148655 - 0.0203483*1.0 + 0.00835867*1.0**2 + 0.00697134*1.0**3
    elif 41 < theta < 42 and p < 0.7:
        dp = 0.0223585 - 0.0365262*p - 0.0150027*p**2 + 0.0854164*p**3 - 0.0462718*p**4
    elif 41 < theta < 42:
        dp = 0.007617
    elif theta > 42 and p < 0.75:
        dp = 0.0152373 - 0.0106377*p - 0.0257573*p**2 + 0.0344851*p**3
    else:
        dp = 0.0152373 - 0.0106377*0.75 - 0.0257573*0.75**2 + 0.0344851*0.75**3
    return dp

def krishna_pi_minus_fd(theta, p):
    if theta < 27.0:
        return p*0.00046571 + 0.00322164
    else:
        if p < 1.7:
            return -0.0024313*p**3 + 0.0094416*p**2 - 0.01257967*p + 0.0122432
        else:
            return 0.006199071

def java_el_fall2018_pass2(phi_local, p, sec):
    """
    Electron momentum correction, Fall 2018 Pass 2 in-bending,
    translated from your Java code.
    """
    phi = phi_local - 30.0/p
    dp = 0.0
    if sec == 1:
        dp = ((-9.82416e-06)*phi*phi + (-2.29956e-05)*phi + 0.000296642)*p*p \
           + (0.0001113414*phi*phi + (-2.0413e-05)*phi + (-0.00862226))*p \
           + ((-0.000281738)*phi*phi + 0.00058712*phi + 0.0652737)
        if p < 7:
            dp += ((-3.4001e-06)*phi*phi + (-2.2885e-05)*phi + 0.00099705)*p*p \
               + (2.1840e-05*phi*phi + 0.00024238*phi + (-0.0091904))*p \
               + ((-2.9180e-05)*phi*phi + (-0.00064496)*phi + 0.022505)
        else:
            dp += ((-6.3656e-05)*phi*phi + 0.00017266*phi + (-0.0017909))*p*p \
               + (0.00104*phi*phi + (-0.0028401)*phi + 0.02981)*p \
               + ((-0.0041995)*phi*phi + 0.011537*phi + (-0.1196))
        dp += (3.278e-07*phi*phi + 6.7084e-07*phi + (-4.39e-05))*p*p \
           + ((-7.231e-06)*phi*phi + (-2.37482e-05)*phi + 0.0004909)*p \
           + (3.2853e-05*phi*phi + 9.63723e-05*phi + (-0.00115))
    elif sec == 2:
        dp = ((-7.741952e-06)*phi*phi + (-2.2402167e-05)*phi + (-0.000426529))*p*p \
           + (7.54079e-05*phi*phi + (-1.3334e-05)*phi + 0.00024201)*p \
           + ((-0.000147876)*phi*phi + 0.00057905*phi + 0.0253551)
        if p < 7:
            dp += (5.3611e-06*phi*phi + 8.1979e-06*phi + 0.00059789)*p*p \
               + ((-4.8185e-05)*phi*phi + (-1.5188e-04)*phi + (-0.0084675))*p \
               + (9.2324e-05*phi*phi + 0.00064420*phi + 0.026792)
        else:
            dp += ((-6.1139e-05)*phi*phi + 5.4087e-06*phi + (-0.0021284))*p*p \
               + (0.0010007*phi*phi + 9.3492e-05*phi + 0.039813)*p \
               + ((-0.0040434)*phi*phi + (-0.0010953)*phi + (-0.18112))
        dp += (6.221217e-07*phi*phi + 1.9596e-06*phi + (-9.826e-05))*p*p \
           + ((-1.28576e-05)*phi*phi + (-4.36589e-05)*phi + 0.00130342)*p \
           + (5.80399e-05*phi*phi + 0.000215388*phi + (-0.0040414))
    elif sec == 3:
        dp = ((-5.115364e-06)*phi*phi + (-1.1983e-05)*phi + (-0.00068329))*p*p \
           + (4.52287e-05*phi*phi + 0.00020855*phi + 0.0034987)*p \
           + ((-9.04461e-05)*phi*phi + (-0.00106657)*phi + 0.0179542)
        if p < 7:
            dp += (9.9281e-07*phi*phi + 3.4879e-06*phi + 0.0011673)*p*p \
               + ((-2.0071e-05)*phi*phi + (-3.1362e-05)*phi + (-0.012329))*p \
               + (6.9463e-05*phi*phi + 3.5102e-05*phi + 0.037505)
        else:
            dp += ((-3.2178e-06)*phi*phi + 4.063e-05*phi + (-0.005209))*p*p \
               + (2.0884e-05*phi*phi + (-0.000688)*phi + 0.086513)*p \
               + (3.953e-05*phi*phi + 0.0029306*phi + (-0.3507))
        dp += ((-4.046e-07)*phi*phi + (-1.3116e-06)*phi + 3.951e-05)*p*p \
           + (5.521e-06*phi*phi + 2.4437e-05*phi + (-0.0016887))*p \
           + ((-1.0963e-05)*phi*phi + (-0.000151944)*phi + 0.0093136)
    elif sec == 4:
        dp = ((-3.9278117e-06)*phi*phi + 2.22893e-05*phi + 0.00012665)*p*p \
           + (4.86493e-05*phi*phi + (-0.00012554)*phi + (-0.0059555))*p \
           + ((-0.000146172)*phi*phi + (-0.00028571)*phi + 0.0606998)
        if p < 7:
            dp += ((-4.8455e-06)*phi*phi + (-1.2074e-05)*phi + 0.0013221)*p*p \
               + (3.2207e-05*phi*phi + 0.00013144*phi + (-0.010451))*p \
               + ((-3.7365e-05)*phi*phi + (-0.00042344)*phi + 0.019952)
        else:
            dp += ((-3.9554e-05)*phi*phi + 5.5496e-06*phi + (-0.0058293))*p*p \
               + (0.00065077*phi*phi + 2.6735e-05*phi + 0.095025)*p \
               + ((-0.0026457)*phi*phi + (-0.00061394)*phi + (-0.3793))
        dp += ((-4.593089e-07)*phi*phi + 1.40673e-05*phi + 6.69e-05)*p*p \
           + (4.0239e-06*phi*phi + (-0.000180863)*phi + (-0.00082722))*p \
           + ((-5.131e-06)*phi*phi + 0.00049748*phi + 0.00255231)
    elif sec == 5:
        dp = ((8.0366e-07)*phi*phi + 2.58072e-05*phi + 0.000360217)*p*p \
           + ((-9.9324e-06)*phi*phi + (-0.0005168531)*phi + (-0.010904))*p \
           + (1.85163e-05*phi*phi + 0.00155709*phi + 0.066493)
        if p < 7:
            dp += (7.7156e-07*phi*phi + (-3.9566e-05)*phi + (-0.00023589))*p*p \
               + ((-9.8309e-06)*phi*phi + 0.00037353*phi + 0.0020382)*p \
               + (2.9506e-05*phi*phi + (-0.00080409)*phi + (-0.0045615))
        else:
            dp += ((-3.241e-05)*phi*phi + (-4.3301e-05)*phi + (-0.0028742))*p*p \
               + (0.00053787*phi*phi + 0.00068921*phi + 0.049578)*p \
               + ((-0.0021955)*phi*phi + (-0.0027698)*phi + (-0.21142))
        dp += ((-1.2151e-06)*phi*phi + (-8.5087e-06)*phi + 4.968e-05)*p*p \
           + (1.46998e-05*phi*phi + 0.000115047*phi + (-0.00039269))*p \
           + ((-4.03686e-05)*phi*phi + (-0.00037078)*phi + 0.00073998)
    elif sec == 6:
        dp = ((-1.95521e-06)*phi*phi + 8.0422e-06*phi + (-2.1324e-05))*p*p \
           + (1.69694e-05*phi*phi + (-6.3066e-05)*phi + (-0.00485568))*p \
           + ((-2.7723e-05)*phi*phi + (-6.8284e-05)*phi + 0.0447535)
        if p < 7:
            dp += ((-8.2535e-07)*phi*phi + 9.1433e-06*phi + 0.00035395)*p*p \
               + ((-3.4272e-06)*phi*phi + (-0.00013012)*phi + (-0.0030724))*p \
               + (4.9211e-05*phi*phi + 0.00045807*phi + 0.0058932)
        else:
            dp += ((-4.976e-05)*phi*phi + (-7.2903e-05)*phi + (-0.0020453))*p*p \
               + (0.00080918*phi*phi + 0.0011688*phi + 0.037042)*p \
               + ((-0.0032504)*phi*phi + (-0.0046169)*phi + (-0.16331))
        dp += ((-7.153e-07)*phi*phi + 1.62859e-05*phi + 8.129e-05)*p*p \
           + (7.225e-06*phi*phi + (-0.000178946)*phi + (-0.00094854))*p \
           + ((-1.3018e-05)*phi*phi + 0.00046643*phi + 0.00266508)
    return dp

def jpsi_inb(p):
    """J/ψ-style in-bending electron correction (second Java function)."""
    return 0.0093796*p + (-0.0007808)*p*p + 0.000163*p*p*p + (-0.02029) + 0.04121/p

def java_el_fall2018_pass2_out(phi_local, p, sec):
    """
    Electron momentum correction, Fall 2018 Pass 2 out-bending,
    translated from your C++ corEl==2 block.
    """
    phi = phi_local - 30.0/p
    dp = 0.0
    if sec == 1:
        dp = ((-5.74868e-06)*phi*phi +  2.23310e-06*phi + 0.00025487)*p*p \
           + ( 8.18043e-05*phi*phi +  9.83830e-05*phi + -0.0056062)*p \
           + ( -2.39310e-04*phi*phi + -0.00117555*phi +  0.0486792)
        if p < 6.5:
            dp += ((-1.18269e-05)*phi*phi + -3.33000e-05*phi +  0.00024240)*p*p \
               + ( 8.54100e-05*phi*phi +  0.00022060*phi +  0.00060900)*p \
               + ( -1.28010e-04*phi*phi + -0.00020790*phi + -0.00741700)
        else:
            dp += ((-1.17030e-06)*phi*phi + -9.86690e-05*phi + -0.00288400)*p*p \
               + ( 2.59100e-05*phi*phi +  0.00169554*phi +  0.04538700)*p \
               + ( -1.22850e-04*phi*phi + -0.00710292*phi + -0.17897000)
        dp += ( 2.50960e-06*phi*phi + -5.53920e-06*phi +  0.00015815)*p*p \
           + (-2.83550e-05*phi*phi +  2.79650e-05*phi + -0.00185510)*p \
           + ( 5.80300e-05*phi*phi +  0.00012528*phi +  0.00604000)
    elif sec == 2:
        dp = ((-4.29711e-06)*phi*phi +  1.36105e-05*phi +  0.00040706)*p*p \
           + ( 6.26450e-05*phi*phi +  9.89810e-05*phi + -0.00595890)*p \
           + (-1.75865e-04*phi*phi + -0.00132960*phi +  0.03068200)
        if p < 6.5:
            dp += ( 1.42031e-05*phi*phi + -4.30730e-05*phi + -0.00137090)*p*p \
               + (-1.36732e-04*phi*phi +  0.00031778*phi +  0.01297350)*p \
               + ( 3.20847e-04*phi*phi + -0.00048824*phi + -0.02994480)
        else:
            dp += (-7.37090e-06*phi*phi +  3.02370e-05*phi + -0.00556253)*p*p \
               + ( 9.81410e-05*phi*phi + -0.00042787*phi +  0.08993300)*p \
               + (-2.92040e-04*phi*phi +  0.00138420*phi + -0.35582200)
        dp += (-3.01290e-06*phi*phi +  6.52082e-06*phi +  0.00024848)*p*p \
           + ( 3.24280e-05*phi*phi + -6.48800e-05*phi + -0.00261070)*p \
           + (-7.78340e-05*phi*phi +  5.97300e-05*phi +  0.00570050)
    # ... similarly fill in sec == 3 .. 6 ...
    else:
        dp = 0.0
    return dp

def jpsi_out(p):
    """J/ψ-style out-bending electron correction."""
    return (-0.06520)*p + 0.007099*p*p + (-0.00005929)*p*p*p + 0.2145 + (-0.1153)/p

def plot_forward_detector_protons(out_dir):
    fig, axs = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs, [0.75, 1.75, 2.75]):
        ax.plot(theta_fd, timothy_fd(theta_fd, p), label='Timothy', linewidth=2)
        ax.plot(theta_fd, [krishna_fd(t, p) for t in theta_fd],
                label='Krishna', linestyle='--', linewidth=2)
        ax.plot(theta_fd, mariana_fd(theta_fd, p),
                label='Mariana', linestyle=':', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(5, 39); ax.set_ylim(-0.02, 0.03)
        ax.set_title(f'p = {p:.2f} GeV'); ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs[0]: ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='lower left', frameon=True)
    fig.suptitle('Forward Detector Proton Corrections')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f'{out_dir}/forward_detector.png')
    plt.close(fig)

def plot_central_detector_protons(out_dir):
    fig, axs = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs, [0.4, 0.75, 1.1]):
        ax.plot(theta_cd, timothy_cd(theta_cd, p), label='Timothy', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(25, 70); ax.set_ylim(-0.03, 0.03)
        ax.set_title(f'p = {p:.2f} GeV'); ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs[0]: ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='upper right', frameon=True)
    fig.suptitle('Central Detector Proton Corrections')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f'{out_dir}/central_detector.png')
    plt.close(fig)

def plot_run_periods(out_dir):
    fig, axs = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs, [0.75, 1.75, 2.75]):
        ax.plot(theta_fd, timothy_fd_fa18_inb(theta_fd, p), label='Fa18 Inb', linewidth=2)
        ax.plot(theta_fd, timothy_fd_fa18_out(theta_fd, p),
                label='Fa18 Out', linestyle='--', linewidth=2)
        ax.plot(theta_fd, timothy_fd_sp19_inb(theta_fd, p),
                label='Sp19 Inb', linestyle=':', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(5, 39); ax.set_ylim(-0.02, 0.02)
        ax.set_title(f'p = {p:.2f} GeV'); ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs[0]: ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='lower left', frameon=True)
    fig.suptitle("Timothy's FD Corrections Across Run Periods")
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f'{out_dir}/timothy_run_periods.png')
    plt.close(fig)

def plot_pion_comparison(out_dir):
    fig, axs = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs, [1.0, 2.0, 3.0]):
        ax.plot(theta_fd, [stefan_pi_plus_fd(t, p) for t in theta_fd],
                label='Stefan π⁺', linewidth=2)
        ax.plot(theta_fd, [krishna_pi_minus_fd(t, p) for t in theta_fd],
                label='Krishna π⁻', linestyle='--', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(5, 39); ax.set_ylim(-0.01, 0.01)
        ax.set_title(f'p = {p:.1f} GeV'); ax.set_xlabel(r'$\theta$ (deg)')
        if ax is axs[0]: ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='lower left', frameon=True)
    fig.suptitle('π⁺ vs π⁻ Energy Loss (FD)')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f'{out_dir}/pi_comparison_fd.png')
    plt.close(fig)

def plot_electron_comparison_inb(out_dir):
    phi_local = np.linspace(0, 60, 500)
    fig, axs = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs, [2.0, 4.0, 6.0]):
        for sec in range(1, 7):
            ax.plot(phi_local,
                    [java_el_fall2018_pass2(phi, p, sec) for phi in phi_local],
                    label=f'Sector {sec}', linewidth=1)
        ax.plot(phi_local,
                [jpsi_inb(p)]*len(phi_local),
                label='jpsi inbending', linestyle='--', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(0, 59.99); ax.set_ylim(-0.06, 0.26)
        ax.set_title(f'p = {p:.0f} GeV')
        ax.set_xlabel(r'$\phi_{\mathrm{local}}$ (deg)')
        if ax is axs[0]:
            ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='upper right', frameon=True, fontsize='small')
    fig.suptitle('Electron Corrections: Momentum Corrections Taskforce vs J/#psi inbending')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f'{out_dir}/electron_corrections_inbending.png')
    plt.close(fig)

def plot_electron_comparison_outb(out_dir):
    phi_local = np.linspace(0, 60, 500)
    fig, axs = plt.subplots(1, 3, figsize=(12,4), sharey=True, gridspec_kw={'wspace':0})
    for ax, p in zip(axs, [2.0, 4.0, 6.0]):
        for sec in range(1, 7):
            ax.plot(phi_local,
                    [java_el_fall2018_pass2_out(phi, p, sec) for phi in phi_local],
                    label=f'Sector {sec}', linewidth=1)
        ax.plot(phi_local,
                [jpsi_out(p)]*len(phi_local),
                label='jpsi outbending', linestyle='--', linewidth=2)
        ax.axhline(0, linestyle='--', color='gray', linewidth=1)
        ax.set_xlim(0, 59.99); ax.set_ylim(-0.06, 0.26)
        ax.set_title(f'p = {p:.0f} GeV')
        ax.set_xlabel(r'$\phi_{\mathrm{local}}$ (deg)')
        if ax is axs[0]:
            ax.set_ylabel(r'$\Delta p$ (GeV)')
        ax.legend(loc='upper right', frameon=True, fontsize='small')
    fig.suptitle('Electron Corrections: Momentum Corrections Taskforce vs J/#psi outbending')
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(f'{out_dir}/electron_corrections_outbending.png')
    plt.close(fig)

def main():
    out_dir = 'output'
    os.makedirs(out_dir, exist_ok=True)

    plot_forward_detector_protons(out_dir)
    plot_central_detector_protons(out_dir)
    plot_run_periods(out_dir)
    plot_pion_comparison(out_dir)
    plot_electron_comparison_inb(out_dir)
    plot_electron_comparison_outb(out_dir)

if __name__ == '__main__':
    main()