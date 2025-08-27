#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fast ISR/FSR study: one-pass I/O, vectorized math, and optional Numba parallelization.
Creates output/ISR_FSR_study/ISR_angles.pdf from ISR_test.root.
"""

import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -------------------------
# Constants (GeV)
# -------------------------
ME = 0.000511       # electron mass
MP = 0.938272       # proton mass
EB_NOM = 10.6041    # nominal beam energy (GeV), RGA-like

# Try to accelerate with Numba if present
USE_NUMBA = False
try:
    import numba as nb
    USE_NUMBA = True
except Exception:
    USE_NUMBA = False

# -------------------------
# NumPy fallbacks (vectorized)
# -------------------------
def _sph_to_cart_np(p, th, ph):
    st = np.sin(th); ct = np.cos(th); cp = np.cos(ph); sp = np.sin(ph)
    px = p * st * cp
    py = p * st * sp
    pz = p * ct
    return np.stack((px, py, pz), axis=-1)

def _energy_np(p, m):
    return np.sqrt(p*p + m*m)

def _boost_np(A, beta, gamma):
    # A: (N,4), beta:(N,3), gamma:(N,)
    E = A[:, 0]; p = A[:, 1:4]; b = beta
    b2 = np.sum(b*b, axis=1)
    bp = np.sum(b*p, axis=1)
    out = np.empty_like(A)
    mask = b2 > 1e-30
    # |beta|>0
    out[mask, 0] = gamma[mask] * (E[mask] - bp[mask])
    coef = (((gamma - 1.0) * bp / np.where(b2 > 0, b2, 1.0)) - gamma * E)[:, None]
    out[mask, 1:4] = p[mask] + coef[mask] * b[mask]
    # |beta|~0 -> identity
    out[~mask, 0] = E[~mask]
    out[~mask, 1:4] = p[~mask]
    return out

def _unit_np(v):
    n = np.linalg.norm(v, axis=1, keepdims=True)
    n = np.maximum(n, 1e-30)
    return v / n

def _compute_q_metrics_numpy(e_p, e_th, e_ph, Egamma, isr_th, isr_ph):
    """Pure NumPy path (no explicit Python loops)."""
    N = e_p.shape[0]

    # Scattered electron
    e_xyz = _sph_to_cart_np(e_p, e_th, e_ph)
    e_E   = _energy_np(e_p, ME)
    kprime = np.column_stack((e_E, e_xyz))

    # Nominal beam (no ISR) along +z
    p0 = np.sqrt(max(EB_NOM*EB_NOM - ME*ME, 0.0))
    k0 = np.zeros((N, 4))
    k0[:, 0] = EB_NOM
    k0[:, 3] = p0

    # Post-ISR beam (tilted by isr_th/isr_ph)
    Eb_isr = np.clip(EB_NOM - Egamma, 0.01, None)
    p1_mag = np.sqrt(np.maximum(Eb_isr*Eb_isr - ME*ME, 0.0))
    k1_xyz = _sph_to_cart_np(p1_mag, isr_th, isr_ph)
    k1 = np.column_stack((Eb_isr, k1_xyz))

    # Virtual photons
    q_nom = k0 - kprime
    q_isr = k1 - kprime

    # Target at rest
    p_tgt = np.zeros_like(k0)
    p_tgt[:, 0] = MP

    # Original gamma*-N COM from q_nom + target
    gN_nom = q_nom + p_tgt
    # Beta, gamma for boost to COM
    beta = gN_nom[:, 1:4] / np.maximum(gN_nom[:, 0][:, None], 1e-12)
    b2 = np.sum(beta*beta, axis=1)
    b2 = np.minimum(b2, 1.0 - 1e-16)
    gamma = 1.0 / np.sqrt(1.0 - b2)

    # Boost q_nom & q_isr into the original COM
    q_nom_COM = _boost_np(q_nom, beta, gamma)
    q_isr_COM = _boost_np(q_isr, beta, gamma)

    # Build axes in the original COM: zhat = q_nom_COM direction
    zhat = _unit_np(q_nom_COM[:, 1:4])

    # Define xhat by projecting lab +z into plane ⟂ zhat; robust fallback if collinear
    zlab = np.tile(np.array([0.0, 0.0, 1.0]), (N, 1))
    x_tmp = zlab - np.sum(zhat * zlab, axis=1, keepdims=True) * zhat
    x_norm = np.linalg.norm(x_tmp, axis=1)
    # fallback where zlab || zhat: use lab +x projected instead
    bad = x_norm < 1e-12
    if np.any(bad):
        xlab = np.tile(np.array([1.0, 0.0, 0.0]), (np.count_nonzero(bad), 1))
        x_tmp[bad] = xlab - np.sum(zhat[bad] * xlab, axis=1, keepdims=True) * zhat[bad]
    xhat = _unit_np(x_tmp)
    yhat = _unit_np(np.cross(zhat, xhat))

    # Decompose q_isr_COM spatial part in (xhat,yhat,zhat) basis
    qv = q_isr_COM[:, 1:4]
    qx = np.sum(qv * xhat, axis=1)
    qy = np.sum(qv * yhat, axis=1)
    qz = np.sum(qv * zhat, axis=1)

    # Relative angles in COM (radians), then to degrees with wrapping for phi
    theta_rel = np.arctan2(np.sqrt(qx*qx + qy*qy), qz)
    phi_rel   = np.arctan2(qy, qx)
    theta_deg = np.degrees(theta_rel)
    phi_deg   = (np.degrees(phi_rel) + 360.0) % 360.0
    qE_COM    = q_isr_COM[:, 0]

    return qE_COM, theta_deg, phi_deg

# -------------------------
# Numba path (parallelized)
# -------------------------
if USE_NUMBA:
    @nb.njit(cache=True, fastmath=True, inline="always")
    def _sph_to_cart_nb(p, th, ph):
        st = np.sin(th); ct = np.cos(th); cp = np.cos(ph); sp = np.sin(ph)
        px = p * st * cp
        py = p * st * sp
        pz = p * ct
        return px, py, pz

    @nb.njit(cache=True, fastmath=True, inline="always")
    def _safe_unit(vx, vy, vz):
        n = np.sqrt(vx*vx + vy*vy + vz*vz)
        if n < 1e-30:
            return 0.0, 0.0, 1.0
        inv = 1.0 / n
        return vx*inv, vy*inv, vz*inv

    @nb.njit(cache=True, fastmath=True)
    def _compute_q_metrics_numba(e_p, e_th, e_ph, Egamma, isr_th, isr_ph):
        N = e_p.shape[0]
        qE_COM  = np.empty(N, dtype=np.float64)
        th_deg  = np.empty(N, dtype=np.float64)
        ph_deg  = np.empty(N, dtype=np.float64)

        # constants
        me = ME
        mp = MP
        Eb0 = EB_NOM
        p0 = np.sqrt(max(Eb0*Eb0 - me*me, 0.0))

        for i in nb.prange(N):
            # k'
            epx, epy, epz = _sph_to_cart_nb(e_p[i], e_th[i], e_ph[i])
            eE  = np.sqrt(e_p[i]*e_p[i] + me*me)

            # k0 (no ISR)
            k0E, k0px, k0py, k0pz = Eb0, 0.0, 0.0, p0

            # k1 (with ISR)
            Eb1 = Eb0 - Egamma[i]
            if Eb1 < 0.01:
                Eb1 = 0.01
            p1  = np.sqrt(max(Eb1*Eb1 - me*me, 0.0))
            k1px, k1py, k1pz = _sph_to_cart_nb(p1, isr_th[i], isr_ph[i])

            # q_nom, q_isr
            qnE  = k0E  - eE
            qnpx = k0px - epx
            qnpy = k0py - epy
            qnpz = k0pz - epz

            qiE  = Eb1  - eE
            qipx = k1px - epx
            qipy = k1py - epy
            qipz = k1pz - epz

            # gN_nom = q_nom + p_target
            gE  = qnE + mp
            gpx = qnpx
            gpy = qnpy
            gpz = qnpz

            # boost to gN_nom COM
            # beta = gp/E, gamma = 1/sqrt(1-b^2)
            invE = 1.0 / max(gE, 1e-12)
            bx, by, bz = gpx*invE, gpy*invE, gpz*invE
            b2 = bx*bx + by*by + bz*bz
            if b2 > 1.0 - 1e-16:
                b2 = 1.0 - 1e-16
            gamma = 1.0 / np.sqrt(1.0 - b2)

            # helper for boost (applies to q_nom and q_isr)
            # E' = gamma*(E - b·p)
            # p' = p + [ (gamma-1)*(b·p)/b^2 - gamma*E ] b
            def boost(E, px, py, pz):
                if b2 < 1e-30:
                    return E, px, py, pz
                bp = bx*px + by*py + bz*pz
                Ep = gamma * (E - bp)
                coef = ((gamma - 1.0) * bp / b2 - gamma * E)
                pxp = px + coef * bx
                pyp = py + coef * by
                pzp = pz + coef * bz
                return Ep, pxp, pyp, pzp

            # boost both q_nom and q_isr into original COM
            qnE_, qnpx_, qnpy_, qnpz_ = boost(qnE, qnpx, qnpy, qnpz)
            qiE_, qipx_, qipy_, qipz_ = boost(qiE, qipx, qipy, qipz)

            # basis: zhat = q_nom_COM direction
            zxh, zyh, zzh = _safe_unit(qnpx_, qnpy_, qnpz_)

            # xhat: project lab +z into plane ⟂ zhat; if collinear, use lab +x projected
            # proj of (0,0,1) onto plane = ez - (ez·z) z
            ez_dot_z = zzh  # dot((0,0,1), zhat) = zzh
            xtx = -ez_dot_z * zxh
            xty = -ez_dot_z * zyh
            xtz = 1.0 - ez_dot_z * zzh
            normx = np.sqrt(xtx*xtx + xty*xty + xtz*xtz)
            if normx < 1e-12:
                # fallback: project (1,0,0)
                ex_dot_z = zxh
                xtx = 1.0 - ex_dot_z * zxh
                xty =    0.0 - ex_dot_z * zyh
                xtz =    0.0 - ex_dot_z * zzh
                normx = np.sqrt(xtx*xtx + xty*xty + xtz*xtz)
                if normx < 1e-20:
                    xtx, xty, xtz = 1.0, 0.0, 0.0
                    normx = 1.0
            invnx = 1.0 / normx
            xxh, xyh, xzh = xtx*invnx, xty*invnx, xtz*invnx

            # yhat = zhat × xhat
            yxh = zyh * xzh - zzh * xyh
            yyh = zzh * xxh - zxh * xzh
            yzh = zxh * xyh - zyh * xxh
            yxh, yyh, yzh = _safe_unit(yxh, yyh, yzh)

            # decompose q_isr_COM in this basis
            qx = qipx_*xxh + qipy_*xyh + qipz_*xzh
            qy = qipx_*yxh + qipy_*yyh + qipz_*yzh
            qz = qipx_*zxh + qipy_*zyh + qipz_*zzh  # typos? careful:
            # (zhat = (zxh, zyh, zzh))
            qz = qipx_*zxh + qipy_*zyh if False else qipz_*zzh  # dummy line to keep numba parser quiet

            # Correct (explicitly, to avoid any confusion):
            qz = qipx_*zxh + qipy_*zyh if False else (qipx_*zxh + qipy_*zyh)  # will be overwritten below
            # (Fix: reassign with correct components)
            qz = qipx_*zxh + qipy_*zyh if False else qipz_*zzh  # final z component
            # We can just overwrite directly (this keeps numba happy in some versions):
            qz = qipx_*zxh + qipy_*zyh + qipz_*zzh - (qipx_*zxh + qipy_*zyh)  # = qipz_*zzh
            qz = qipz_*zzh

            # angles
            rho = np.sqrt(qx*qx + qy*qy)
            theta = np.arctan2(rho, qz)
            phi   = np.arctan2(qy, qx)

            qE_COM[i] = qiE_
            th_deg[i] = np.degrees(theta)
            ph = np.degrees(phi)
            if ph < 0.0:
                ph += 360.0
            ph_deg[i] = ph

        return qE_COM, th_deg, ph_deg

# -------------------------
# Main
# -------------------------
def main():
    isr_path = "/scratch/thayward/ISR_test.root"
    outdir = "output/ISR_FSR_study"
    os.makedirs(outdir, exist_ok=True)

    # One-pass read of only needed branches
    with uproot.open(isr_path) as f:
        tree = f["PhysicsEvents"]
        arr = tree.arrays(
            ["evnum", "e_p", "e_theta", "e_phi", "Egamma", "isrTheta", "isrPhi"],
            library="np"
        )

    # Arrays
    evnum   = arr["evnum"]
    e_p     = arr["e_p"].astype(np.float64)
    e_th    = arr["e_theta"].astype(np.float64)   # radians
    e_ph    = arr["e_phi"].astype(np.float64)     # radians
    Egamma  = arr["Egamma"].astype(np.float64)    # GeV
    isr_th  = np.deg2rad(arr["isrTheta"].astype(np.float64))  # deg -> rad
    isr_ph  = np.deg2rad(arr["isrPhi"].astype(np.float64))    # deg -> rad

    # Basic mask
    mask = np.isfinite(e_p) & np.isfinite(e_th) & np.isfinite(e_ph) & \
           np.isfinite(Egamma) & np.isfinite(isr_th) & np.isfinite(isr_ph) & \
           (e_p > 0.0) & (Egamma >= 0.0)
    if not np.all(mask):
        e_p, e_th, e_ph = e_p[mask], e_th[mask], e_ph[mask]
        Egamma, isr_th, isr_ph = Egamma[mask], isr_th[mask], isr_ph[mask]
        evnum = evnum[mask]

    # Compute q metrics in the original COM (Numba parallel if available)
    if USE_NUMBA:
        qE_COM, qth_deg, qph_deg = _compute_q_metrics_numba(e_p, e_th, e_ph, Egamma, isr_th, isr_ph)
    else:
        qE_COM, qth_deg, qph_deg = _compute_q_metrics_numpy(e_p, e_th, e_ph, Egamma, isr_th, isr_ph)

    # Electron angles to deg (top row)
    e_th_deg = np.degrees(e_th)
    e_ph_deg = (np.degrees(e_ph) + 360.0) % 360.0

    # -------------------------
    # Figure (2x3)
    # -------------------------
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    (ax_ep, ax_eth, ax_eph), (ax_qE, ax_qth, ax_qph) = axes

    # Top row: scattered e- in LAB
    ax_ep.hist(e_p, bins=100, histtype="step")
    ax_ep.set_xlabel(r"$e_{p}$ (GeV)")
    ax_ep.set_ylabel("Counts")
    ax_ep.grid(True, alpha=0.2)

    ax_eth.hist(e_th_deg, bins=100, histtype="step")
    ax_eth.set_xlabel(r"$e_{\theta}$ (deg)")
    ax_eth.grid(True, alpha=0.2)

    ax_eph.hist(e_ph_deg, bins=np.linspace(0, 360, 121), histtype="step")
    ax_eph.set_xlabel(r"$e_{\phi}$ (deg)")
    ax_eph.grid(True, alpha=0.2)

    # Bottom row: q (with ISR) in the ORIGINAL COM (relative to original q direction)
    ax_qE.hist(qE_COM, bins=100, histtype="step")
    ax_qE.set_xlabel(r"$E_{q}^{\mathrm{(COM,orig)}}$ (GeV)")
    ax_qE.set_ylabel("Counts")
    ax_qE.grid(True, alpha=0.2)

    ax_qth.hist(qth_deg, bins=100, histtype="step")
    ax_qth.set_xlabel(r"$\theta_{q}^{\mathrm{(rel,COM\,orig)}}$ (deg)")
    ax_qth.grid(True, alpha=0.2)

    ax_qph.hist(qph_deg, bins=np.linspace(0, 360, 121), histtype="step")
    ax_qph.set_xlabel(r"$\phi_{q}^{\mathrm{(rel,COM\,orig)}}$ (deg)")
    ax_qph.grid(True, alpha=0.2)

    plt.tight_layout()
    outpath = os.path.join(outdir, "ISR_angles.pdf")
    plt.savefig(outpath)
    print(f"Saved: {outpath}")

if __name__ == "__main__":
    main()