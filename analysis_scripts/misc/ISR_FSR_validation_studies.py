#!/usr/bin/env python3
import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

# -------------------------
# Constants (GeV, radians)
# -------------------------
ME = 0.000511      # electron mass
MP = 0.938272      # proton mass
EB_NOM = 10.6041   # nominal beam energy (GeV)

# -------------------------
# Small helpers (vectorized)
# -------------------------
def sph_to_cart(p, theta, phi):
    """
    p in GeV, theta, phi in radians. Returns (px, py, pz) with z = beam axis.
    Shapes broadcast along the last dimension.
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
    Given 4-vector(s) P = (E, px, py, pz), return beta (N,3) and gamma (N,).
    """
    E = P[..., 0]
    pvec = P[..., 1:4]
    beta = pvec / np.expand_dims(np.maximum(E, 1e-12), axis=-1)
    b2 = np.sum(beta*beta, axis=-1)
    b2 = np.minimum(b2, 1.0 - 1e-16)  # keep numeric safety
    gamma = 1.0 / np.sqrt(1.0 - b2)
    return beta, gamma

def lorentz_boost(A, beta, gamma):
    """
    Boost 4-vector(s) A = (E, px, py, pz) by velocity -beta (i.e., into the rest frame of P if beta,gamma came from P).
    Vectorized. Formula:
      p' = p + [ (gamma-1)*(b·p)/b^2 - gamma*E ] * b
      E' = gamma*(E - b·p)
    """
    E = A[..., 0]
    p = A[..., 1:4]
    b = beta
    b2 = np.sum(b*b, axis=-1)
    bp = np.sum(b * p, axis=-1)

    out = np.empty_like(A)
    mask = b2 > 1e-30

    # where |beta|>0
    Eprime = np.empty_like(E)
    Eprime[mask] = gamma[mask] * (E[mask] - bp[mask])
    coef = ((gamma - 1.0) * bp / np.where(b2 > 0, b2, 1.0) - gamma * E)[..., None]
    pprime = np.empty_like(p)
    pprime[mask] = p[mask] + coef[mask] * b[mask]

    # where |beta| ~ 0, identity
    Eprime[~mask] = E[~mask]
    pprime[~mask] = p[~mask]

    out[..., 0] = Eprime
    out[..., 1:4] = pprime
    return out

def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    n = np.maximum(n, 1e-30)
    return v / n

def angle_theta_phi(v, xhat, yhat, zhat):
    """
    Spherical angles of v in the orthonormal basis (xhat,yhat,zhat):
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
born_path     = "/scratch/thayward/Born_test.root"
isr_path      = "/scratch/thayward/ISR_test.root"
isr_fsr_path  = "/scratch/thayward/ISR_FSR_test.root"

outdir = "output/ISR_FSR_study"
os.makedirs(outdir, exist_ok=True)

# Only ISR file is needed for this figure
with uproot.open(isr_path) as f:
    tree = f["PhysicsEvents"]
    arr = tree.arrays(
        ["evnum",
         "e_p", "e_theta", "e_phi",      # scattered electron (LAB)
         "Egamma", "isrTheta", "isrPhi"  # ISR info (angles stored in degrees)
        ],
        library="np"
    )

# -------------------------
# Prepare inputs (vectorized)
# -------------------------
evnum  = arr["evnum"].astype(np.int64)

# scattered electron in LAB (your Java saved these in radians)
e_p    = arr["e_p"].astype(np.float64)
e_th   = arr["e_theta"].astype(np.float64)
e_ph   = arr["e_phi"].astype(np.float64)

# ISR information
Egamma = arr["Egamma"].astype(np.float64)               # GeV
isr_th = np.deg2rad(arr["isrTheta"].astype(np.float64)) # radians
isr_ph = np.deg2rad(arr["isrPhi"].astype(np.float64))   # radians

# Basic quality mask
mask = np.isfinite(e_p) & np.isfinite(e_th) & np.isfinite(e_ph) \
       & np.isfinite(Egamma) & np.isfinite(isr_th) & np.isfinite(isr_ph) \
       & (e_p > 0.0) & (Egamma >= 0.0)
e_p, e_th, e_ph = e_p[mask], e_th[mask], e_ph[mask]
Egamma, isr_th, isr_ph = Egamma[mask], isr_th[mask], isr_ph[mask]
evnum = evnum[mask]

N = e_p.shape[0]

# -------------------------
# 4-vectors
# -------------------------
# Scattered electron k' (LAB)
kprime_E     = energy_from_p_mass(e_p, ME)
kprime_pxyz  = sph_to_cart(e_p, e_th, e_ph)
kprime       = np.column_stack((kprime_E, kprime_pxyz))  # (N,4)

# Nominal (no ISR) beam k0 along +z
p0 = np.sqrt(max(EB_NOM*EB_NOM - ME*ME, 0.0))
k0 = np.zeros((N, 4))
k0[:, 0] = EB_NOM
k0[:, 3] = p0

# Post-ISR beam k1: energy EB_NOM - Egamma, direction (isr_th, isr_ph)
Eb_isr   = np.clip(EB_NOM - Egamma, 0.01, None)
p1       = np.sqrt(np.maximum(Eb_isr*Eb_isr - ME*ME, 0.0))
k1_pxyz  = sph_to_cart(p1, isr_th, isr_ph)
k1       = np.column_stack((Eb_isr, k1_pxyz))           # (N,4)

# --- Top row uses the *incident* beam after ISR ---
beam_p       = p1
beam_th_deg  = np.rad2deg(isr_th)
beam_ph_deg  = (np.rad2deg(isr_ph) + 360.0) % 360.0

# Virtual photons q = k - k'
q_nom = k0 - kprime
q_isr = k1 - kprime

# Compute Q^2 = -q^2 with metric (+,-,-,-)
def Q2_from_q(q):
    q2 = q[:, 0]*q[:, 0] - np.sum(q[:, 1:4]*q[:, 1:4], axis=-1)
    return np.maximum(0.0, -q2)  # guard tiny negatives from roundoff

Q2_nom = Q2_from_q(q_nom)
Q2_isr = Q2_from_q(q_isr)
dQ2    = Q2_isr - Q2_nom  # Delta Q2 = (with ISR) - (no ISR)

# -------------------------
# Angles of q^{ISR} in the ORIGINAL gamma*-N COM frame
# -------------------------
# Target proton at rest (LAB)
p_target = np.zeros_like(k0)
p_target[:, 0] = MP

# Define original COM by (q_nom + p)
gN_nom         = q_nom + p_target
beta0, gamma0  = boost_vector(gN_nom)         # boost to COM of (q_nom+p)

# Boost both q_nom and q_isr into that same COM
q_nom_COM = lorentz_boost(q_nom, beta0, gamma0)  # reference for axes
q_isr_COM = lorentz_boost(q_isr, beta0, gamma0)  # to be measured

# Build orthonormal basis from the original q direction
z_hat = unit(q_nom_COM[:, 1:4])
z_lab = np.tile(np.array([0.0, 0.0, 1.0]), (N, 1))      # +z in LAB
x_tmp = z_lab - np.sum(z_hat * z_lab, axis=-1, keepdims=True) * z_hat
x_hat = unit(x_tmp)
y_hat = unit(np.cross(z_hat, x_hat))

# Spherical angles of q_isr in that basis
q_isr_COM_p         = q_isr_COM[:, 1:4]
theta_rel, phi_rel  = angle_theta_phi(q_isr_COM_p, x_hat, y_hat, z_hat)
theta_rel_deg       = np.rad2deg(theta_rel)
phi_rel_deg         = (np.rad2deg(phi_rel) + 360.0) % 360.0

# -------------------------
# Plot: 2x3 canvas
# -------------------------
fig, axes = plt.subplots(2, 3, figsize=(12, 7))
(ax_ep, ax_eth, ax_eph), (ax_dQ2, ax_qth, ax_qph) = axes

# Top row: *incident beam after ISR* (not the scattered electron)
# -- make y-axis log and relabel x to k_p (GeV)
n_ep, _, _ = ax_ep.hist(beam_p, bins=100, histtype="step")
ax_ep.set_xlabel(r"$k_p$ (GeV)")
ax_ep.set_ylabel("Counts")
ax_ep.set_yscale("log")
pos = n_ep[n_ep > 0]
if pos.size > 0:
    ax_ep.set_ylim(bottom=max(1.0, 0.8 * pos.min()))

ax_eth.hist(beam_th_deg, bins=100, histtype="step")
ax_eth.set_xlabel(r"$k_{\theta}$ (deg)")

ax_eph.hist(beam_ph_deg, bins=np.linspace(0, 360, 121), histtype="step")
ax_eph.set_xlabel(r"$k_{\phi}$ (deg)")

# Bottom-left: Delta Q^2 = Q^2_ISR - Q^2_nom
ax_dQ2.hist(dQ2, bins=120, histtype="step")
ax_dQ2.set_xlabel(r"$\Delta Q^2 \equiv Q^2_{\rm ISR}-Q^2_{\rm nom}$ (GeV$^2$)")
ax_dQ2.set_ylabel("Counts")

# Bottom-middle / right: angles of q^{ISR} in the ORIGINAL COM basis
ax_qth.hist(theta_rel_deg, bins=100, histtype="step")
ax_qth.set_xlabel(r"$\theta_{q}^{\mathrm{(rel,COM\,orig)}}$ (deg)")

ax_qph.hist(phi_rel_deg, bins=np.linspace(0, 360, 121), histtype="step")
ax_qph.set_xlabel(r"$\phi_{q}^{\mathrm{(rel,COM\,orig)}}$ (deg)")

for ax in axes.flat:
    ax.grid(True, alpha=0.2)
#endfor

plt.tight_layout()
outpath = os.path.join(outdir, "ISR_angles.pdf")
plt.savefig(outpath)
print(f"Saved: {outpath}")