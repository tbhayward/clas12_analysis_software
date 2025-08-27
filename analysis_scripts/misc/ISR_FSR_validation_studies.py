#!/usr/bin/env python3
import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

# -------------------------
# Constants (GeV, radians)
# -------------------------
ME = 0.000511   # electron mass
MP = 0.938272   # proton mass
EB_NOM = 10.6041  # nominal beam energy (GeV)

# -------------------------
# Small helpers
# -------------------------
def sph_to_cart(p, theta, phi):
    """
    p in GeV, theta, phi in radians. Returns (px, py, pz) with z = beam axis.
    """
    st = np.sin(theta)
    ct = np.cos(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)
    px = p * st * cp
    py = p * st * sp
    pz = p * ct
    return np.stack((px, py, pz), axis=-1)

def energy_from_p_mass(p, m):
    return np.sqrt(p*p + m*m)

def boost_vector(P):
    """
    Given a 4-vector P = (E, px, py, pz), return beta (3,) and gamma (scalar).
    """
    E = P[..., 0]
    pvec = P[..., 1:4]
    beta = pvec / np.expand_dims(np.maximum(E, 1e-12), axis=-1)
    b2 = np.sum(beta*beta, axis=-1)
    b2 = np.minimum(b2, 1.0 - 1e-16)
    gamma = 1.0 / np.sqrt(1.0 - b2)
    return beta, gamma

def lorentz_boost(A, beta, gamma):
    """
    Boost 4-vector(s) A = (E, px, py, pz) by velocity -beta (i.e., into the rest frame of P if beta,gamma came from P).
    Vectorized over first dimension.
    Formula: p' = p + [ (gamma-1)*(b·p)/b^2 - gamma*E ] * b
             E' = gamma*(E - b·p)
    """
    # shapes: A:(N,4), beta:(N,3), gamma:(N,)
    E = A[..., 0]
    p = A[..., 1:4]
    b = beta
    b2 = np.sum(b*b, axis=-1)
    # handle b2 ~ 0 safely
    mask = b2 > 1e-30
    bp = np.sum(b * p, axis=-1)

    # allocate outputs
    Eprime = np.empty_like(E)
    pprime = np.empty_like(p)

    # where |beta|>0
    Eprime[mask] = gamma[mask]*(E[mask] - bp[mask])
    coef = ((gamma - 1.0) * bp / np.where(b2 > 0, b2, 1.0) - gamma * E)[..., None]
    pprime[mask] = p[mask] + coef[mask] * b[mask]

    # where |beta| ~ 0, it's identity
    Eprime[~mask] = E[~mask]
    pprime[~mask] = p[~mask]

    out = np.empty_like(A)
    out[..., 0] = Eprime
    out[..., 1:4] = pprime
    return out

def unit(v):
    """Return unit vector(s) with safe handling of zero-norm rows."""
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    n = np.maximum(n, 1e-30)
    return v / n

def angle_theta_phi(v, xhat, yhat, zhat):
    """
    Given a vector v and an orthonormal basis (xhat, yhat, zhat), compute spherical angles:
      theta = arctan2(sqrt(vx^2+vy^2), vz)
      phi   = atan2(vy, vx)
    Returns (theta, phi) in radians.
    """
    vx = np.sum(v * xhat, axis=-1)
    vy = np.sum(v * yhat, axis=-1)
    vz = np.sum(v * zhat, axis=-1)
    theta = np.arctan2(np.sqrt(vx*vx + vy*vy), vz)
    phi = np.arctan2(vy, vx)
    return theta, phi

# -------------------------
# I/O
# -------------------------
born_path = "/scratch/thayward/Born_test.root"
isr_path = "/scratch/thayward/ISR_test.root"
isr_fsr_path = "/scratch/thayward/ISR_FSR_test.root"

outdir = "output/ISR_FSR_study"
os.makedirs(outdir, exist_ok=True)

# Only ISR_test.root is needed for the first figure
with uproot.open(isr_path) as f:
    tree = f["PhysicsEvents"]

    # Load just what we need now
    arr = tree.arrays([
        "evnum",
        "e_p", "e_theta", "e_phi",
        "Egamma", "isrTheta", "isrPhi"
    ], library="np")

evnum   = arr["evnum"]

# Electron kinematics (scattered e- in LAB)
e_p     = arr["e_p"]                # GeV/c
e_th    = arr["e_theta"]            # radians (from your Java)
e_ph    = arr["e_phi"]              # radians (0..2π)

# ISR info
Egamma  = arr["Egamma"]             # GeV
# saved in degrees -> convert to radians for math
isr_th  = np.deg2rad(arr["isrTheta"])
isr_ph  = np.deg2rad(arr["isrPhi"])

# Basic quality masks (avoid nonsense)
mask = np.isfinite(e_p) & np.isfinite(e_th) & np.isfinite(e_ph) & np.isfinite(Egamma) \
       & np.isfinite(isr_th) & np.isfinite(isr_ph) \
       & (e_p > 0.0) & (Egamma >= 0.0)
if not np.all(mask):
    e_p, e_th, e_ph = e_p[mask], e_th[mask], e_ph[mask]
    Egamma, isr_th, isr_ph = Egamma[mask], isr_th[mask], isr_ph[mask]
    evnum = evnum[mask]

# ---------------------------------------------
# Build 4-vectors needed for both rows
# ---------------------------------------------

# Scattered electron 4-vector k' (LAB)
e_pxpyz = sph_to_cart(e_p, e_th, e_ph)                       # (N,3)
e_E     = energy_from_p_mass(e_p, ME)
kprime  = np.column_stack((e_E, e_pxpyz))                    # (N,4)

# Nominal (no-ISR) beam 4-vector k0 along +z
p0 = np.sqrt(max(EB_NOM*EB_NOM - ME*ME, 0.0))
k0 = np.zeros((kprime.shape[0], 4))
k0[:, 0] = EB_NOM
k0[:, 3] = p0  # pz

# Post-ISR beam 4-vector k1 tilted by (isr_th, isr_ph), with Eb_isr = EB_NOM - Egamma
Eb_isr = np.clip(EB_NOM - Egamma, 0.01, None)
p1 = np.sqrt(np.maximum(Eb_isr*Eb_isr - ME*ME, 0.0))
k1_pxpyz = sph_to_cart(p1, isr_th, isr_ph)
k1 = np.column_stack((Eb_isr, k1_pxpyz))

# Virtual photons
q_nom = k0 - kprime        # no-ISR virtual photon (reference)
q_isr = k1 - kprime        # virtual photon with ISR

# Target proton at rest (LAB)
p_target = np.zeros_like(k0)
p_target[:, 0] = MP

# ---------------------------------------------
# Define the ORIGINAL gamma*-N COM frame (no ISR)
# ---------------------------------------------
gN_nom = q_nom + p_target              # total system 4-p for original case
beta0, gamma0 = boost_vector(gN_nom)   # boost to bring gN_nom to rest

# Boost q_nom and q_isr into that original COM frame
q_nom_COM = lorentz_boost(q_nom, beta0, gamma0)
q_isr_COM = lorentz_boost(q_isr, beta0, gamma0)

# Build the "original" axis system in that COM:
#   z-hat = unit(q_nom_COM spatial)
#   x-hat = projection of lab +z (0,0,1) onto plane ⟂ z-hat
#   y-hat = z-hat × x-hat
z_hat = unit(q_nom_COM[:, 1:4])

z_lab = np.tile(np.array([0.0, 0.0, 1.0]), (z_hat.shape[0], 1))
x_tmp = z_lab - np.sum(z_hat * z_lab, axis=-1, keepdims=True) * z_hat
x_hat = unit(x_tmp)
y_hat = unit(np.cross(z_hat, x_hat))

# Measure q_isr in that basis
q_isr_COM_p = q_isr_COM[:, 1:4]
q_isr_COM_E = q_isr_COM[:, 0]
theta_rel, phi_rel = angle_theta_phi(q_isr_COM_p, x_hat, y_hat, z_hat)  # radians

# Convert everything to degrees for plotting
e_th_deg  = np.rad2deg(e_th)
e_ph_deg  = np.rad2deg(e_ph)
theta_rel_deg = np.rad2deg(theta_rel)
phi_rel_deg   = (np.rad2deg(phi_rel) + 360.0) % 360.0

# -------------------------
# Make the 2x3 figure
# -------------------------
fig, axes = plt.subplots(2, 3, figsize=(12, 7))
(ax_ep, ax_eth, ax_eph), (ax_qE, ax_qth, ax_qph) = axes

# Top row: electron lab momentum and angles (histograms)
ax_ep.hist(e_p, bins=100, histtype="step")
ax_ep.set_xlabel(r"$e_{p}$ (GeV)")
ax_ep.set_ylabel("Counts")

ax_eth.hist(e_th_deg, bins=100, histtype="step")
ax_eth.set_xlabel(r"$e_{\theta}$ (deg)")

ax_eph.hist(e_ph_deg, bins=100, histtype="step")
ax_eph.set_xlabel(r"$e_{\phi}$ (deg)")

# Bottom row: q (with ISR) measured in the ORIGINAL COM frame
ax_qE.hist(q_isr_COM_E, bins=100, histtype="step")
ax_qE.set_xlabel(r"$E_{q}^{\mathrm{(COM,orig)}}$ (GeV)")
ax_qE.set_ylabel("Counts")

ax_qth.hist(theta_rel_deg, bins=100, histtype="step")
ax_qth.set_xlabel(r"$\theta_{q}^{\mathrm{(rel,COM\,orig)}}$ (deg)")

ax_qph.hist(phi_rel_deg, bins=np.linspace(0,360,121), histtype="step")
ax_qph.set_xlabel(r"$\phi_{q}^{\mathrm{(rel,COM\,orig)}}$ (deg)")

for ax in axes.flat:
    ax.grid(True, alpha=0.2)

plt.tight_layout()
outpath = os.path.join(outdir, "ISR_angles.pdf")
plt.savefig(outpath)
print(f"Saved: {outpath}")