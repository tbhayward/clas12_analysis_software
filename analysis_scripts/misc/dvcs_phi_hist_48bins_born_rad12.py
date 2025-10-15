import os
import math
import argparse
import glob
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROTON_MASS  = 0.938272

def unit(v):
    n = np.linalg.norm(v)
    if n == 0.0:
        return v
    #end if
    return v / n
#end def

def fourvec(E, px, py, pz):
    return np.array([E, px, py, pz], dtype=float)
#end def

def lorentz_boost(p4, beta_vec):
    # Boost 4-vector p4 = (E, px, py, pz) by 3-velocity beta vector.
    bx, by, bz = beta_vec
    b2 = bx*bx + by*by + bz*bz
    if b2 <= 0.0:
        return p4
    #end if
    gamma = 1.0 / math.sqrt(1.0 - b2)
    bp = bx*p4[1] + by*p4[2] + bz*p4[3]
    gamma2 = (gamma - 1.0) / b2

    E_prime  = gamma * (p4[0] - bp)
    px_prime = p4[1] + (-gamma * bx * p4[0]) + gamma2 * bp * bx
    py_prime = p4[2] + (-gamma * by * p4[0]) + gamma2 * bp * by
    pz_prime = p4[3] + (-gamma * bz * p4[0]) + gamma2 * bp * bz
    return np.array([E_prime, px_prime, py_prime, pz_prime], dtype=float)
#end def

def parse_lund_events_streaming(path, max_events=None):
    """
    Streaming generator over events from a LUND text file.

    dict:
      - beam_E 
      - target_id 
      - weight 
      - particles: list of dict with keys: pid, E, px, py, pz, mass
    """
    count = 0
    with open(path, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            #end if
            header = header.strip()
            if len(header) == 0:
                continue
            #end if

            toks = header.split()
            if len(toks) < 10:
                continue
            #end if

            try:
                Npart = int(float(toks[0]))
            except Exception:
                continue
            #end if

            # column 7 beam energy, 8 interacted nucleon id, 10 weight
            try:
                beam_E    = float(toks[6])
                target_id = int(float(toks[7]))
                weight    = float(toks[9])
            except Exception:
                beam_E    = float(toks[6]) if len(toks) > 6 else float("nan")
                target_id = 2212
                weight    = 1.0
            #end if

            particles = []
            for _ in range(Npart):
                pline = f.readline()
                if not pline:
                    break
                #end if
                pline = pline.strip()
                ptok = pline.split()
                if len(ptok) < 14:
                    continue
                #end if

                try:
                    pid = int(float(ptok[3]))
                    px  = float(ptok[6])
                    py  = float(ptok[7])
                    pz  = float(ptok[8])
                    E   = float(ptok[9])
                    m   = float(ptok[10])
                except Exception:
                    continue
                #end if

                particles.append({"pid": pid, "E": E, "px": px, "py": py, "pz": pz, "mass": m})
            #end for

            if len(particles) != Npart:
                continue
            #end if

            yield {
                "beam_E": beam_E,
                "target_id": target_id,
                "weight": weight,
                "particles": particles,
            }

            count += 1
            if max_events is not None and count >= max_events:
                break
            #end if
        #end while
    #end with
#end def

# Photon selectors 

def select_gamma_born_only_one(particles):
    """
    For Born: return only photon
    """
    photons = [p for p in particles if p["pid"] == 22]
    if len(photons) == 1:
        return photons[0]
    #end if
    return None
#end def

def select_gamma_first(particles):
    """
    For Rad-first: return the first photon
    """
    for p in particles:
        if p["pid"] == 22:
            return p
        #end if
    #end for
    return None
#end def

def select_gamma_second_or_first(particles):
    """
    For Rad-second-or-first: second photon if exists, otherwise the first.
    """
    seen_first = None
    for p in particles:
        if p["pid"] == 22:
            if seen_first is None:
                seen_first = p
            else:
                return p  # second
            #end if
        #end if
    #end for
    return seen_first
#end def

def trento_phi_in_gammaN_cm(beam_E, target_id, particles, gamma_selector):
    """
    Compute Trento phi for the selected photon in the gamma*-N CM frame.
    gamma_selector: function(particles) -> photon dict or None
    Returns (phi, ok_flag). Phi in [0, 2*pi).
    """
    # Incoming e- along +z
    k = fourvec(beam_E, 0.0, 0.0, beam_E)

    eles = [p for p in particles if p["pid"] == 11]
    if len(eles) != 1:
        return (float("nan"), False)
    #end if
    eprime = eles[0]
    kprime = fourvec(eprime["E"], eprime["px"], eprime["py"], eprime["pz"])

    gamma = gamma_selector(particles)
    if gamma is None:
        return (float("nan"), False)
    #end if
    r = fourvec(gamma["E"], gamma["px"], gamma["py"], gamma["pz"])

    MN = PROTON_MASS
    P  = fourvec(MN, 0.0, 0.0, 0.0)

    q = k - kprime
    if q[0] <= 0.0:
        return (float("nan"), False)
    #end if

    # Boost to gamma*-N CM
    W = q + P
    beta_vec  = W[1:4] / W[0]
    k_cm      = lorentz_boost(k,      beta_vec)
    kprime_cm = lorentz_boost(kprime, beta_vec)
    q_cm      = lorentz_boost(q,      beta_vec)
    r_cm      = lorentz_boost(r,      beta_vec)

    k_vec      = k_cm[1:4]
    kprime_vec = kprime_cm[1:4]
    q_vec      = q_cm[1:4]
    r_vec      = r_cm[1:4]

    z_q   = unit(q_vec)
    n_lep = unit(np.cross(k_vec, kprime_vec))
    n_had = unit(np.cross(q_vec, r_vec))

    # if np.linalg.norm(n_lep) == 0.0 or np.linalg.norm(n_had) == 0.0 or np.linalg.norm(z_q) == 0.0:
    #     return (float("nan"), False)
    # #end if

    y_val = np.dot(np.cross(n_lep, n_had), z_q)
    x_val = np.dot(n_lep, n_had)
    phi   = math.atan2(y_val, x_val)
    if phi < 0.0:
        phi += 2.0 * math.pi
    #endif
    return (phi, True)
#end def

# ------------------------------- Directory -> histogram with TOTAL event cap -------------------------------

def list_matching_files(dir_path, pattern):
    if not os.path.isdir(dir_path):
        raise RuntimeError("Directory not found: {}".format(dir_path))
    #endif
    files = glob.glob(os.path.join(dir_path, pattern))
    files = sorted([f for f in files if os.path.isfile(f)])
    return files
#enddef

def hist_from_directory_totalcap(dir_path, nbins, use_weights, pattern,
                                 max_files, max_total_events, gamma_selector):
    """
    Returns (hist_counts, edges, n_events_used)
    Reads across files but stops once n_events_used reaches max_total_events.
    """
    files = list_matching_files(dir_path, pattern)
    if len(files) == 0:
        raise RuntimeError("No files matched in {} with pattern {}".format(dir_path, pattern))
    #endif
    if max_files is not None:
        files = files[:max_files]
    #endif

    edges  = np.linspace(0.0, 2.0 * math.pi, nbins + 1)
    counts = np.zeros(nbins, dtype=float)
    twopi  = 2.0 * math.pi

    n_events_used = 0
    for path in files:
        if max_total_events is not None and n_events_used >= max_total_events:
            break
        #endif
        # Stream and stop when we hit the total cap
        for evt in parse_lund_events_streaming(path, max_events=None):
            if max_total_events is not None and n_events_used >= max_total_events:
                break
            #endif
            phi, ok = trento_phi_in_gammaN_cm(evt["beam_E"], evt["target_id"], evt["particles"], gamma_selector)
            if not ok or not np.isfinite(phi):
                continue
            #endif
            idx = int(nbins * (phi / twopi))
            if idx == nbins:
                idx = nbins - 1
            #endif
            incr = evt["weight"] if use_weights else 1.0
            counts[idx] += incr
            n_events_used += 1
        #endfor
    #endfor

    return counts, edges, n_events_used
#enddef

# ------------------------------- Plotting -------------------------------

def plot_overlay_counts_log(edges, h_born, h_rad_first, h_rad_second_or_first, out_png, title):
    # Replace zeros with NaN to avoid log(0) issues
    def safe_for_log(y):
        return np.where(y > 0.0, y, np.nan)
    #enddef

    hb   = safe_for_log(h_born)
    hr1  = safe_for_log(h_rad_first)
    hr21 = safe_for_log(h_rad_second_or_first)

    centers = 0.5 * (edges[:-1] + edges[1:])

    plt.figure(figsize=(9, 5.5))
    plt.yscale("log")
    plt.plot(centers, hb,   marker="o", linestyle="-", label="Born")
    plt.plot(centers, hr1,  marker="s", linestyle="--", label="Rad (first gamma)")
    plt.plot(centers, hr21, marker="^", linestyle=":", label="Rad (second gamma or first)")
    plt.xlabel("phi_trento (rad)", fontsize=13)
    plt.ylabel("raw counts (log scale)", fontsize=13)
    plt.title(title, fontsize=14)
    plt.xlim(0.0, 2.0 * math.pi)
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
#enddef

# ------------------------------- Main -------------------------------

def main():
    ap = argparse.ArgumentParser(description="Make 100-bin RAW-COUNT phi histograms with a total event cap per dataset (Born and Rad).")
    ap.add_argument("--born-dir", required=True, help="Directory containing Born LUND files")
    ap.add_argument("--rad-dir",  required=True, help="Directory containing Radiative LUND files")
    ap.add_argument("--outdir", default="output/dvcs_phi_hist_100bins_raw", help="Output directory")
    ap.add_argument("--pattern-born", default="*.dat", help="Glob pattern for Born files (default *.dat)")
    ap.add_argument("--pattern-rad",  default="*.dat", help="Glob pattern for Rad files (default *.dat)")
    ap.add_argument("--use-weights", action="store_true", help="Use LUND event weights from header")
    ap.add_argument("--max-files-born", type=int, default=None, help="Optional limit on number of Born files")
    ap.add_argument("--max-files-rad",  type=int, default=None, help="Optional limit on number of Rad files")
    ap.add_argument("--max-total-events", type=int, default=100000, help="Total event cap per dataset (default 100000)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    nbins = 100

    print("======================================")
    print("[INFO] Born dir           :", args.born_dir)
    print("[INFO] Rad  dir           :", args.rad_dir)
    print("[INFO] Patterns (born,rad): {}, {}".format(args.pattern_born, args.pattern_rad))
    print("[INFO] nbins              :", nbins)
    print("[INFO] weights            :", "LUND header weights" if args.use_weights else "unit weights")
    print("[INFO] max files (born,rad): {}, {}".format(args.max_files_born, args.max_files_rad))
    print("[INFO] max total events    :", args.max_total_events)
    print("======================================")

    # 1) Born photon: only events with exactly one photon
    h_born, edges, n_born = hist_from_directory_totalcap(
        args.born_dir, nbins, args.use_weights, args.pattern_born,
        args.max_files_born, args.max_total_events, gamma_selector=select_gamma_born_only_one
    )

    # 2) First photon in the rad file
    h_rad_first, _, n_rad_first = hist_from_directory_totalcap(
        args.rad_dir, nbins, args.use_weights, args.pattern_rad,
        args.max_files_rad, args.max_total_events, gamma_selector=select_gamma_first
    )

    # 3) Second photon if exists, else first photon (rad)
    h_rad_second_or_first, _, n_rad_2or1 = hist_from_directory_totalcap(
        args.rad_dir, nbins, args.use_weights, args.pattern_rad,
        args.max_files_rad, args.max_total_events, gamma_selector=select_gamma_second_or_first
    )

    out_png = os.path.join(args.outdir, "phi_hist_100bins_raw_born_rad_first_second_or_first_log.png")
    title = "DVCS phi_Trento (gamma*-N CM), 100 bins — RAW counts (log scale)"
    plot_overlay_counts_log(edges, h_born, h_rad_first, h_rad_second_or_first, out_png, title)

    print("[INFO] Events used (capped): Born={}, Rad-first={}, Rad-(2nd or 1st)={}".format(
        n_born, n_rad_first, n_rad_2or1))
    print("[INFO] Sum counts          : Born={:.0f}, Rad-first={:.0f}, Rad-(2nd or 1st)={:.0f}".format(
        np.sum(h_born), np.sum(h_rad_first), np.sum(h_rad_second_or_first)))
    print("[SAVED]", out_png)
#enddef

if __name__ == "__main__":
    main()
#endif