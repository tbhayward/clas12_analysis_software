#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make a 2x2 canvas from two LUND .dat files (Born vs Radiative) for exclusive DVCS:
  [1,1]  electron momentum e_p (GeV)
  [1,2]  Q^2 (GeV^2)
  [2,1]  W (GeV)
  [2,2]  -t (GeV^2)  (t is the Mandelstam variable)

Assumptions:
- Beam: 10.6 GeV electron along +z on a stationary proton target.
- LUND columns (1-based):
    header: at least 10 values; the 1st is N (particles per event)
    particle lines (per event):
        1: index
        2: lifetime (ns)  [unused]
        3: type (only type==1 is propagated)  [unused here]
        4: PDG id
        5-6: parent/first-daughter [unused]
        7-9: px, py, pz  (GeV)
        10: E (GeV)      [UD; may or may not be present/accurate]
        11: mass (GeV)   [UD; may be present]
        12-14: vx, vy, vz (cm) [unused]

Usage:
    python lund_rad_plots.py born.dat radiative.dat

Output:
    /u/home/thayward/rad_test.pdf
"""

import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt

# Use plain ASCII minus in tick labels (avoids Unicode minus).
plt.rcParams["axes.unicode_minus"] = False

E_BEAM = 10.6  # (GeV)
M_E = 0.000511 # (GeV)
M_P = 0.9382720813  # proton mass (GeV)

# Minimal PDG mass table for common particles in DVCS final state
PDG_MASS = {
    11: M_E,        # e-
    -11: M_E,       # e+
    22: 0.0,        # gamma
    2212: M_P,      # proton
    2112: 0.939565, # neutron (not expected in ep->epgamma, but harmless)
}

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")
# endif

def energy_from_pm(px, py, pz, mass):
    p2 = px*px + py*py + pz*pz
    m2 = max(mass*mass, 0.0)
    return math.sqrt(max(p2 + m2, 0.0))
# endif

def four_sq(e, px, py, pz):
    # Minkowski metric (+,-,-,-)
    return e*e - px*px - py*py - pz*pz
# endif

def lund_event_stream(path):
    """
    Generator over events. Yields a list of particle dicts with keys:
      pdg, px, py, pz, E, m
    Robust to extra spaces/blank lines.
    """
    with open(path, "r") as f:
        line_iter = iter(f)
        for line in line_iter:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) < 1:
                continue
            # Heuristic: header line must have >= 1 value and first is N (int)
            # LUND specs say header has >= 10 values; be tolerant of spacing.
            if len(parts) >= 1:
                # Attempt to parse N and require len(parts) >= 1 (we'll accept any header >= 1 and < 100)
                try:
                    N = int(float(parts[0]))
                except Exception:
                    continue
                # Read N particle lines
                particles = []
                read_count = 0
                while read_count < N:
                    try:
                        pline = next(line_iter)
                    except StopIteration:
                        break
                    ps = pline.strip()
                    if not ps:
                        continue
                    cols = ps.split()
                    # Require at least up to pz (col 9)
                    if len(cols) < 9:
                        continue
                    # 1-based -> 0-based indices
                    pdg = int(float(cols[3])) if len(cols) >= 4 else 0
                    px  = safe_float(cols[6]) if len(cols) >= 7 else float("nan")
                    py  = safe_float(cols[7]) if len(cols) >= 8 else float("nan")
                    pz  = safe_float(cols[8]) if len(cols) >= 9 else float("nan")
                    E   = safe_float(cols[9]) if len(cols) >= 10 else float("nan")
                    m   = safe_float(cols[10]) if len(cols) >= 11 else float("nan")

                    # Fallbacks for mass / energy if missing or invalid
                    if not math.isfinite(m):
                        m = PDG_MASS.get(pdg, 0.0)
                    if not math.isfinite(E) or E <= 0.0:
                        E = energy_from_pm(px, py, pz, m)

                    particles.append({"pdg": pdg, "px": px, "py": py, "pz": pz, "E": E, "m": m})
                    read_count += 1
                # endfor

                if particles:
                    yield particles
                # endif
            # endif
        # endfor
    # with
# enddef

def event_kinematics(particles, ebeam=E_BEAM):
    """
    Extract e', p' and gamma' from particles, compute:
      - e_p (|p_e'|)
      - Q2  using q = k - k'
      - W   via (p + q)^2 with target at rest
      - t   preferring hadronic definition (p' - p)^2, else photon-based (q - q')^2
    Returns tuple (e_p, Q2, W, minus_t) or None if not computable.
    """
    # Find final-state objects (pick the highest-momentum electron if multiple)
    electrons = [p for p in particles if p["pdg"] == 11]
    protons   = [p for p in particles if p["pdg"] == 2212]
    photons   = [p for p in particles if p["pdg"] == 22]

    if not electrons:
        return None
    # endif

    eprime = max(electrons, key=lambda a: a["px"]*a["px"] + a["py"]*a["py"] + a["pz"]*a["pz"])
    ep_px, ep_py, ep_pz, ep_E = eprime["px"], eprime["py"], eprime["pz"], eprime["E"]
    e_p = math.sqrt(ep_px*ep_px + ep_py*ep_py + ep_pz*ep_pz)

    # Beam e- 4-vector
    kE = ebeam
    k  = (kE, 0.0, 0.0, kE)

    # Scattered e- 4-vector
    kp = (ep_E, ep_px, ep_py, ep_pz)

    # q = k - k'
    qE  = k[0] - kp[0]
    qpx = k[1] - kp[1]
    qpy = k[2] - kp[2]
    qpz = k[3] - kp[3]

    q2 = four_sq(qE, qpx, qpy, qpz)
    Q2 = -q2  # positive

    if not (math.isfinite(Q2) and Q2 > 0.0):
        return None
    # endif

    # Target proton at rest
    p_targ = (M_P, 0.0, 0.0, 0.0)

    # W^2 = (p + q)^2
    W2 = four_sq(p_targ[0] + qE, p_targ[1] + qpx, p_targ[2] + qpy, p_targ[3] + qpz)
    W = math.sqrt(max(W2, 0.0))

    # t from recoiling proton if present; else photon-based
    t_val = float("nan")
    if protons:
        pprime = max(protons, key=lambda a: a["px"]*a["px"] + a["py"]*a["py"] + a["pz"]*a["pz"])
        pE, ppx, ppy, ppz = pprime["E"], pprime["px"], pprime["py"], pprime["pz"]
        # (p' - p)^2
        t_val = four_sq(pE - p_targ[0], ppx - p_targ[1], ppy - p_targ[2], ppz - p_targ[3])
    # endif

    if not math.isfinite(t_val) and photons:
        g = max(photon s, key=lambda a: a["px"]*a["px"] + a["py"]*a["py"] + a["pz"]*a["pz"])
        gE, gpx, gpy, gpz = g["E"], g["px"], g["py"], g["pz"]
        # (q - q')^2
        t_val = four_sq(qE - gE, qpx - gpx, qpy - gpy, qpz - gpz)
    # endif

    if not math.isfinite(t_val):
        return None
    # endif

    minus_t = -t_val
    return (e_p, Q2, W, minus_t)
# enddef

def collect_kinematics(path):
    e_p, Q2, W, minus_t = [], [], [], []
    n_ev = 0
    n_ok = 0
    for particles in lund_event_stream(path):
        n_ev += 1
        kin = event_kinematics(particles, ebeam=E_BEAM)
        if kin is None:
            continue
        # endif
        ep, q2, w, mt = kin
        e_p.append(ep)
        Q2.append(q2)
        W.append(w)
        minus_t.append(mt)
        n_ok += 1
    # endfor
    return {
        "e_p": np.array(e_p, dtype=float),
        "Q2": np.array(Q2, dtype=float),
        "W": np.array(W, dtype=float),
        "-t": np.array(minus_t, dtype=float),
        "n_events": n_ev,
        "n_kept": n_ok,
    }
# enddef

def robust_range(a, lo=1.0, hi=99.0, pad=0.05):
    if a.size == 0:
        return (0.0, 1.0)
    # endif
    lo_v = np.nanpercentile(a, lo)
    hi_v = np.nanpercentile(a, hi)
    if not np.isfinite(lo_v) or not np.isfinite(hi_v) or lo_v >= hi_v:
        return (float(np.nanmin(a)), float(np.nanmax(a)))
    # endif
    span = hi_v - lo_v
    return (max(lo_v - pad*span, 0.0), hi_v + pad*span)
# enddef

def main():
    if len(sys.argv) != 3:
        print("Usage: python lund_rad_plots.py <born.dat> <radiative.dat>")
        sys.exit(2)
    # endif
    f_born = sys.argv[1]
    f_rad  = sys.argv[2]

    born = collect_kinematics(f_born)
    rad  = collect_kinematics(f_rad)

    # Combine for consistent axis ranges
    def join(a, b): 
        return np.concatenate([a, b]) if a.size and b.size else (a if a.size else b)
    # enddef

    ep_all = join(born["e_p"], rad["e_p"])
    q2_all = join(born["Q2"],  rad["Q2"])
    w_all  = join(born["W"],   rad["W"])
    mt_all = join(born["-t"],  rad["-t"])

    ep_rng = robust_range(ep_all)
    q2_rng = robust_range(q2_all)
    w_rng  = robust_range(w_all)
    mt_rng = robust_range(mt_all)

    # Figure
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    ax11, ax12, ax21, ax22 = axes[0,0], axes[0,1], axes[1,0], axes[1,1]

    bins = 60

    # Helper to draw one panel
    def draw(ax, data_born, data_rad, rng, xlabel):
        ax.hist(data_born, bins=bins, range=rng, histtype="step", linewidth=1.5, label="Born")
        ax.hist(data_rad,  bins=bins, range=rng, histtype="step", linewidth=1.5, label="Radiative")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("counts")
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.2)
    # enddef

    draw(ax11, born["e_p"], born["e_p"]*0+np.nan, ep_rng, r"$e_{p}$ (GeV)")
    ax11.cla()
    # Redraw with both properly:
    ax11.hist(born["e_p"], bins=bins, range=ep_rng, histtype="step", linewidth=1.5, label="Born")
    ax11.hist(rad["e_p"],  bins=bins, range=ep_rng, histtype="step", linewidth=1.5, label="Radiative")
    ax11.set_xlabel(r"$e_{p}$ (GeV)")
    ax11.set_ylabel("counts")
    ax11.legend(frameon=False)
    ax11.grid(True, alpha=0.2)

    draw(ax12, born["Q2"], rad["Q2"], q2_rng, r"$Q^{2}$ (GeV^{2})")
    draw(ax21, born["W"],  rad["W"],  w_rng,  r"$W$ (GeV)")
    draw(ax22, born["-t"], rad["-t"], mt_rng, r"$-t$ (GeV^{2})")

    fig.tight_layout()

    outpath = "/u/home/thayward/rad_test.pdf"
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath)
    print("Saved:", outpath)
    print("Born kept/total:", born["n_kept"], "/", born["n_events"])
    print("Radiative kept/total:", rad["n_kept"], "/", rad["n_events"])
# enddef

if __name__ == "__main__":
    main()
# endif