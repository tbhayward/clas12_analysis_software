#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# dvcs_phi_ratio_both.py
#
# Build 12-bin phi histograms and norad/rad count ratios from:
#   (A) LUND text files (parsing event-by-event, computing Trento phi)
#   (B) ROOT files with precomputed branch "phi2" in tree "PhysicsEvents"
#
# Outputs:
#   output/dvcs_phi_ratio/
#       phi_counts_lund.csv
#       phi_counts_root.csv
#       phi_ratio_both.png
#
# Usage (defaults to your provided file paths):
#   python3 dvcs_phi_ratio_both.py
#
# or explicitly:
#   python3 dvcs_phi_ratio_both.py \
#     --lund-norad /work/clas12/yijie/Simulation/DVCS/dvcsgen/norad_gen/nor001M1.dat \
#     --lund-rad   /work/clas12/yijie/Simulation/DVCS/dvcsgen/rad_gen/rad001M1.dat \
#     --root-norad /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_10604MeV.root \
#     --root-rad   /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_out_10604MeV.root \
#     --tree PhysicsEvents \
#     --phi-branch phi2 \
#     --root-max-entries 100000
#
import os
import math
import argparse
import numpy as np

# plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# uproot for ROOT I/O
try:
    import uproot
except Exception as e:
    uproot = None
#endif


# ------------------------------- Helpers -------------------------------

def unit(v):
    n = np.linalg.norm(v)
    if n == 0.0:
        return v
    #endif
    return v / n
#enddef

def lorentz_boost(p4, beta_vec):
    """
    Boost a 4-vector p4 = (E, px, py, pz) by 3-velocity beta_vec.
    """
    bx, by, bz = beta_vec
    b2 = bx*bx + by*by + bz*bz
    if b2 <= 0.0:
        return p4
    #endif
    gamma = 1.0 / math.sqrt(1.0 - b2)
    bp = bx*p4[1] + by*p4[2] + bz*p4[3]
    gamma2 = (gamma - 1.0) / b2

    E_prime = gamma * (p4[0] - bp)
    px_prime = p4[1] + (-gamma * bx * p4[0]) + gamma2 * bp * bx
    py_prime = p4[2] + (-gamma * by * p4[0]) + gamma2 * bp * by
    pz_prime = p4[3] + (-gamma * bz * p4[0]) + gamma2 * bp * bz
    return np.array([E_prime, px_prime, py_prime, pz_prime], dtype=float)
#enddef

def fourvec(E, px, py, pz):
    return np.array([E, px, py, pz], dtype=float)
#enddef


# ------------------------------- LUND parsing and phi_Trento -------------------------------

PROTON_MASS_DEFAULT = 0.938272
NEUTRON_MASS_DEFAULT = 0.939565

def parse_lund_events(path):
    """
    Generator over events from a LUND text file.

    Yields dict with:
      - beam_E (float)
      - target_id (int)
      - weight (float)
      - particles: list of dict with keys: pid, E, px, py, pz, mass
    """
    with open(path, "r") as f:
        lines = f.readlines()
    #endif

    i = 0
    nlines = len(lines)
    while i < nlines:
        line = lines[i].strip()
        if len(line) == 0:
            i += 1
            continue
        #endif

        toks = line.split()
        if len(toks) < 10:
            i += 1
            continue
        #endif

        # header column 1 is number of particles
        try:
            Npart = int(float(toks[0]))
        except Exception:
            i += 1
            continue
        #endif

        # column 7 beam energy, 8 interacted nucleon id, 10 weight
        try:
            beam_E = float(toks[6])
            target_id = int(float(toks[7]))
            weight = float(toks[9])
        except Exception:
            beam_E = float(toks[6]) if len(toks) > 6 else float("nan")
            target_id = 2212
            weight = 1.0
        #endif

        particles = []
        i += 1
        for j in range(Npart):
            if i >= nlines:
                break
            #endif
            pline = lines[i].strip()
            ptok = pline.split()
            if len(ptok) < 14:
                i += 1
                continue
            #endif

            # columns: 4 pid, 7 px, 8 py, 9 pz, 10 E, 11 mass
            try:
                pid = int(float(ptok[3]))
                px = float(ptok[6])
                py = float(ptok[7])
                pz = float(ptok[8])
                E  = float(ptok[9])
                m  = float(ptok[10])
            except Exception:
                i += 1
                continue
            #endif
            particles.append({"pid": pid, "E": E, "px": px, "py": py, "pz": pz, "mass": m})
            i += 1
        #endfor

        if len(particles) != Npart:
            continue
        #endif

        yield {
            "beam_E": beam_E,
            "target_id": target_id,
            "weight": weight,
            "particles": particles,
        }
    #endwhile
#enddef

def pick_recoil_mass(particles, target_id):
    if target_id == 2212:
        m_default = PROTON_MASS_DEFAULT
    else:
        m_default = NEUTRON_MASS_DEFAULT
    #endif
    for p in particles:
        if p["pid"] in (2212, 2112) and p["mass"] > 0.0:
            return p["mass"]
        #endif
    #endfor
    return m_default
#enddef

def select_dvcs_photon(particles):
    photons = [p for p in particles if p["pid"] == 22]
    if len(photons) == 0:
        return None
    #endif
    photons.sort(key=lambda p: p["E"], reverse=True)
    return photons[0]
#enddef

def pick_scattered_electron(particles):
    eles = [p for p in particles if p["pid"] == 11]
    if len(eles) == 0:
        return None
    #endif
    eles.sort(key=lambda p: p["E"], reverse=True)
    return eles[0]
#enddef

def trento_phi_in_gammaN_cm(beam_E, target_id, particles):
    """
    Compute Trento phi for the DVCS photon in the gamma*-N CM frame.
    Returns (phi, ok_flag).
    """
    k = fourvec(beam_E, 0.0, 0.0, beam_E)

    eprime = pick_scattered_electron(particles)
    if eprime is None:
        return (float("nan"), False)
    #endif
    kprime = fourvec(eprime["E"], eprime["px"], eprime["py"], eprime["pz"])

    dvcs_photon = select_dvcs_photon(particles)
    if dvcs_photon is None:
        return (float("nan"), False)
    #endif
    r = fourvec(dvcs_photon["E"], dvcs_photon["px"], dvcs_photon["py"], dvcs_photon["pz"])

    MN = pick_recoil_mass(particles, target_id)
    P = fourvec(MN, 0.0, 0.0, 0.0)

    q = k - kprime
    if q[0] <= 0.0:
        return (float("nan"), False)
    #endif

    W = q + P
    beta_vec = W[1:4] / W[0]
    k_cm      = lorentz_boost(k,      beta_vec)
    kprime_cm = lorentz_boost(kprime, beta_vec)
    q_cm      = lorentz_boost(q,      beta_vec)
    r_cm      = lorentz_boost(r,      beta_vec)

    k_vec = k_cm[1:4]
    kprime_vec = kprime_cm[1:4]
    q_vec = q_cm[1:4]
    r_vec = r_cm[1:4]

    z_q = unit(q_vec)
    n_lep = unit(np.cross(k_vec, kprime_vec))
    n_had = unit(np.cross(q_vec, r_vec))

    if np.linalg.norm(n_lep) == 0.0 or np.linalg.norm(n_had) == 0.0 or np.linalg.norm(z_q) == 0.0:
        return (float("nan"), False)
    #endif

    y_val = np.dot(np.cross(n_lep, n_had), z_q)
    x_val = np.dot(n_lep, n_had)
    phi = math.atan2(y_val, x_val)
    if phi < 0.0:
        phi += 2.0 * math.pi
    #endif

    return (phi, True)
#enddef

def counts_from_lund(path, nbins=12, use_weights=False):
    counts = np.zeros(nbins, dtype=float)
    for evt in parse_lund_events(path):
        phi, ok = trento_phi_in_gammaN_cm(evt["beam_E"], evt["target_id"], evt["particles"])
        if not ok or not np.isfinite(phi):
            continue
        #endif
        idx = int(nbins * (phi / (2.0 * math.pi)))
        if idx == nbins:
            idx = nbins - 1
        #endif
        incr = evt["weight"] if use_weights else 1.0
        counts[idx] += incr
    #endfor
    return counts
#enddef


# ------------------------------- ROOT reading -------------------------------

def counts_from_root_phi(root_path, tree_name="PhysicsEvents", phi_branch="phi2", nbins=12, max_entries=100000):
    if uproot is None:
        raise RuntimeError("uproot is not available; please `pip install uproot`.")
    #endif

    with uproot.open(root_path) as f:
        tree = f[tree_name]
        n = min(max_entries, tree.num_entries)
        arr = tree[phi_branch].array(entry_stop=n)
        phi_vals = np.asarray(arr, dtype=float)
    #endwith

    # map to [0, 2*pi)
    twopi = 2.0 * math.pi
    phi_vals = np.mod(phi_vals, twopi)

    edges = np.linspace(0.0, twopi, nbins + 1)
    counts, _ = np.histogram(phi_vals, bins=edges)
    return counts
#enddef


# ------------------------------- Driver -------------------------------

def save_csv(csv_path, centers, counts_norad, counts_rad, ratio):
    with open(csv_path, "w") as out:
        out.write("bin,phi_center_rad,norad_counts,rad_counts,ratio\n")
        for i in range(len(centers)):
            val = ratio[i] if np.isfinite(ratio[i]) else ""
            out.write(f"{i},{centers[i]:.6f},{counts_norad[i]:.0f},{counts_rad[i]:.0f},{val}\n")
        #endfor
    #endif
#enddef

def main():
    ap = argparse.ArgumentParser(description="Compute 12-bin DVCS phi ratios from LUND and ROOT (phi2), plot both on one figure.")
    # LUND inputs
    ap.add_argument("--lund-norad", default="/work/clas12/yijie/Simulation/DVCS/dvcsgen/norad_gen/nor001M1.dat", help="LUND norad file")
    ap.add_argument("--lund-rad",   default="/work/clas12/yijie/Simulation/DVCS/dvcsgen/rad_gen/rad001M1.dat", help="LUND rad file")
    ap.add_argument("--lund-use-weights", action="store_true", help="Use LUND header event weights")
    # ROOT inputs
    ap.add_argument("--root-norad", default="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_10604MeV.root", help="ROOT norad file")
    ap.add_argument("--root-rad",   default="/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_out_10604MeV.root", help="ROOT rad file")
    ap.add_argument("--tree", default="PhysicsEvents", help="Tree name")
    ap.add_argument("--phi-branch", default="phi2", help="Phi branch name in ROOT tree (default: phi2)")
    ap.add_argument("--root-max-entries", type=int, default=100000, help="Max entries from each ROOT file")
    # General
    ap.add_argument("--outdir", default="output/dvcs_phi_ratio", help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    nbins = 12
    edges = np.linspace(0.0, 2.0 * math.pi, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # ---------- LUND ----------
    counts_lund_norad = counts_from_lund(args.lund_norad, nbins=nbins, use_weights=args.lund_use_weights)
    counts_lund_rad   = counts_from_lund(args.lund_rad,   nbins=nbins, use_weights=args.lund_use_weights)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_lund = counts_lund_norad / counts_lund_rad
    #endif
    save_csv(os.path.join(args.outdir, "phi_counts_lund.csv"), centers, counts_lund_norad, counts_lund_rad, ratio_lund)

    # ---------- ROOT (phi2) ----------
    counts_root_norad = counts_from_root_phi(args.root_norad, tree_name=args.tree, phi_branch=args.phi_branch, nbins=nbins, max_entries=args.root_max_entries)
    counts_root_rad   = counts_from_root_phi(args.root_rad,   tree_name=args.tree, phi_branch=args.phi_branch, nbins=nbins, max_entries=args.root_max_entries)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_root = counts_root_norad / counts_root_rad
    #endif
    save_csv(os.path.join(args.outdir, "phi_counts_root.csv"), centers, counts_root_norad, counts_root_rad, ratio_root)

    # ---------- Plot both ----------
    plt.figure(figsize=(8, 5))
    plt.plot(centers, ratio_lund, marker="o", linestyle="-", label="LUND norad / rad")
    plt.plot(centers, ratio_root, marker="s", linestyle="--", label=f"ROOT norad / rad (first {args.root_max_entries})")
    plt.xlabel("phi_trento (rad)", fontsize=13)
    plt.ylabel("count ratio", fontsize=13)
    plt.title("DVCS phi_trento ratio in 12 bins (LUND vs ROOT phi2)", fontsize=14)
    plt.xlim(0.0, 2.0 * math.pi)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    fig_path = os.path.join(args.outdir, "phi_ratio_both.png")
    plt.savefig(fig_path, dpi=200)
    plt.close()

    print("Saved:", os.path.join(args.outdir, "phi_counts_lund.csv"))
    print("Saved:", os.path.join(args.outdir, "phi_counts_root.csv"))
    print("Saved:", fig_path)
#enddef

if __name__ == "__main__":
    main()
#endif