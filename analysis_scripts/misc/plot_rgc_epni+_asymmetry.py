#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

# Toggle for quick debugging (process only first 5 runs per period)
QUICK_RUN = True

# Maximum run number to include
MAX_RUNNUM = 17768

# Path to the CSV of run charges
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output directory and filename
OUTPUT_DIR         = "output/enpi+"
FOUR_PANEL_OUTPUT  = os.path.join(OUTPUT_DIR, "four_panel_Mx2_xB_comparison.pdf")

# File lists: (filepath, label)
NH3_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_NH3_epi+.root", "Su22-NH3"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+.root", "Fa22-NH3"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_NH3_epi+.root", "Sp23-NH3"),
]
C_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_C_epi+.root", "Su22-C"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_C_epi+.root", "Fa22-C"),
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_C_epi+.root", "Sp23-C"),
]
H2_FILES = [
    ("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_sp19_inb_H2_epi+.root", "Sp19-H2"),
]

# -----------------------------------------------------------------------------
# PHYSICAL CONSTANTS
# -----------------------------------------------------------------------------

m_e    = 0.000511   # electron mass (GeV)
m_pi   = 0.13957    # charged pion mass (GeV)
m_p    = 0.938272   # proton mass (GeV)
m_n    = 0.939565   # neutron mass (GeV)

# -----------------------------------------------------------------------------
# DILUTION‐FACTOR DATA FOR BOTTOM‐RIGHT PANEL
# -----------------------------------------------------------------------------

# Su22
enpichi2FitsALUsinphi_Su22 = [
    [0.094266731, 1.184739999, 0.126228827],
    [0.168622194, 0.082358508, 0.012946927],
    [0.254887003, 0.115297657, 0.008854289],
    [0.348247742, 0.144713802, 0.009273995],
    [0.441277333, 0.117431442, 0.012353808],
    [0.535009222, 0.094723130, 0.021626204]
]
x_Su22      = np.array([row[0] for row in enpichi2FitsALUsinphi_Su22])
dil_Su22    = np.array([0.711119, 0.376873, 0.390790, 0.401177, 0.410563, 0.416077])
dil_err_Su22= np.array([0.173083, 0.0175818, 0.00942491, 0.00801711, 0.0108668, 0.0225813])

# Fa22
enpichi2FitsALUsinphi_Fa22 = [
    [0.094825648, 1.221850467, 0.095501230],
    [0.168249360, 0.081525475, 0.009989312],
    [0.255099359, 0.118300134, 0.007065332],
    [0.348334763, 0.140055462, 0.007378878],
    [0.441201661, 0.140928523, 0.009773392],
    [0.534977263, 0.110737425, 0.016852914]
]
x_Fa22      = np.array([row[0] for row in enpichi2FitsALUsinphi_Fa22])
dil_Fa22    = np.array([0.374532, 0.406101, 0.386122, 0.397153, 0.415526, 0.428239])
dil_err_Fa22= np.array([0.117283, 0.00608866, 0.00354573, 0.00301186, 0.00400862, 0.00824607])

# Sp23 (we’ll leave placeholders for future)
x_Sp23      = np.array([])   # to be filled later
dil_Sp23    = np.array([])
dil_err_Sp23= np.array([])

# -----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

def parse_run_charges(csv_path):
    """
    Read clas12_run_info.csv and return a dict { runnum: charge }, skipping runs > MAX_RUNNUM.
    Lines beginning with '#' are ignored.
    """
    run_charges = {}
    if not os.path.exists(csv_path):
        print(f"[ERROR] CSV not found: {csv_path}")
        return run_charges

    with open(csv_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split(',')
            if len(parts) < 2:
                continue
            try:
                runnum = int(parts[0])
                charge = float(parts[1])
            except ValueError:
                continue
            if runnum > MAX_RUNNUM:
                continue
            run_charges[runnum] = charge

    print(f"[DEBUG] Parsed {len(run_charges)} runs ≤ {MAX_RUNNUM} from CSV.")
    return run_charges


def get_beam_energy(runnums):
    """
    Given a NumPy array of run numbers, return a NumPy array of beam energies Eb (in GeV):
      • 6616 ≤ run ≤ 6783 → Eb = 10.1998         (RGA Spring 2019, H₂ data)
      • 16042 ≤ run ≤ 17065 → Eb = 10.5473        (RGC Summer 2022)
      • 17067 ≤ run ≤ 17724 → Eb = 10.5563        (RGC Fall 2022)
      • 17725 ≤ run ≤ 17811 → Eb = 10.5593        (RGC Sp 2023)
    Runs outside these ranges → Eb = 0 (will be filtered out).
    """
    Eb = np.zeros_like(runnums, dtype=float)
    mask_h2 = (runnums >= 6616) & (runnums <= 6783)
    mask1   = (runnums >= 16042) & (runnums <= 17065)
    mask2   = (runnums >= 17067) & (runnums <= 17724)
    mask3   = (runnums >= 17725) & (runnums <= 17811)
    Eb[mask_h2] = 10.1998
    Eb[mask1]   = 10.5473
    Eb[mask2]   = 10.5563
    Eb[mask3]   = 10.5593
    return Eb


def compute_t_array(run_arr, e_p_arr, e_th_arr, e_ph_arr, p_p_arr, p_th_arr, p_ph_arr):
    """
    Compute Mandelstam t event-by-event for ep → e' n π⁺, assuming exclusivity.
    t = (q - p_π)², where q = p_beam - p_e'.

    - run_arr:  array of run numbers
    - e_p_arr:  array of scattered-electron momentum magnitudes (GeV)
    - e_th_arr: array of electron polar angles (radians)
    - e_ph_arr: array of electron azimuthal angles (radians)
    - p_p_arr:  array of pion momentum magnitudes (GeV)
    - p_th_arr: array of pion polar angles (radians)
    - p_ph_arr: array of pion azimuthal angles (radians)

    Returns t_array (in GeV²).
    """
    Eb_arr = get_beam_energy(run_arr)

    # 1) Scattered electron e' 4-vector
    E_e   = np.sqrt(e_p_arr**2 + m_e**2)
    sin_e = np.sin(e_th_arr)
    cos_e = np.cos(e_th_arr)
    ex =  e_p_arr * sin_e * np.cos(e_ph_arr)
    ey =  e_p_arr * sin_e * np.sin(e_ph_arr)
    ez =  e_p_arr * cos_e

    # 2) Pion 4-vector
    E_pi = np.sqrt(p_p_arr**2 + m_pi**2)
    sin_p = np.sin(p_th_arr)
    cos_p = np.cos(p_th_arr)
    px =  p_p_arr * sin_p * np.cos(p_ph_arr)
    py =  p_p_arr * sin_p * np.sin(p_ph_arr)
    pz =  p_p_arr * cos_p

    # 3) Virtual photon q 4-vector: q = p_beam - p_e'
    E_q = Eb_arr - E_e
    qx  = -ex
    qy  = -ey
    qz  = Eb_arr - ez

    # 4) Δ = q - p_π
    dE = E_q - E_pi
    dx = qx - px
    dy = qy - py
    dz = qz - pz

    # 5) t = (ΔE)² - (Δp)²
    t_array = dE**2 - (dx**2 + dy**2 + dz**2)
    return t_array


def process_file(args):
    """
    Worker for parallel histogramming:
    Returns (label, hist_t), where hist_t is the normalized Mx² histogram
    after applying (run ≤ MAX_RUNNUM, Mx² ≤ 2.0, |t| < 1) for that file.
    """
    filepath, label, run_charges, quick = args

    # 1) Open TTree
    tree = uproot.open(filepath)["PhysicsEvents"]

    # 2) Extract arrays
    run_arr  = tree["runnum"].array(library="np").astype(int)
    mx2_arr  = tree["Mx2"].array(library="np")
    e_p_arr  = tree["e_p"].array(library="np")
    e_th_arr = tree["e_theta"].array(library="np")
    e_ph_arr = tree["e_phi"].array(library="np")
    p_p_arr  = tree["p_p"].array(library="np")
    p_th_arr = tree["p_theta"].array(library="np")
    p_ph_arr = tree["p_phi"].array(library="np")

    # 3) Base filter: run ≤ MAX_RUNNUM & Mx² ≤ 2.0
    mask_base = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
    run_base  = run_arr[mask_base]
    mx2_base  = mx2_arr[mask_base]
    e_p_base  = e_p_arr[mask_base]
    e_th_base = e_th_arr[mask_base]
    e_ph_base = e_ph_arr[mask_base]
    p_p_base  = p_p_arr[mask_base]
    p_th_base = p_th_arr[mask_base]
    p_ph_base = p_ph_arr[mask_base]

    # Count how many survive the base filter
    num_base = run_base.size
    print(f"[DEBUG][{label}] base‐filtered events: {num_base}")

    # 4) Compute t and apply |t| < 1 cut
    if num_base > 0:
        t_vals = compute_t_array(
            run_base,
            e_p_base, e_th_base, e_ph_base,
            p_p_base, p_th_base, p_ph_base
        )
        mask_t = np.abs(t_vals) < 1.0
    else:
        mask_t = np.array([], dtype=bool)

    run_t = run_base[mask_t]
    mx2_t = mx2_base[mask_t]

    num_t = run_t.size
    print(f"[DEBUG][{label}] after |t|<1: {num_t} events  (survival {100*num_t/num_base:.1f}% )")

    # Prepare bins
    bins = np.linspace(0.0, 1.5, 101)
    hist_t = np.zeros(len(bins)-1, dtype=float)

    if num_t == 0:
        return (label, hist_t)

    # 5) Select unique runs among those that survive
    if quick:
        seen = []
        for r in run_t:
            if r not in seen:
                seen.append(r)
            if len(seen) >= 5:
                break
        selected_runs = set(seen)
    else:
        selected_runs = set(np.unique(run_t))

    # 6) Total charge Q_t for those runs
    Q_t = 0.0
    for r in selected_runs:
        q = run_charges.get(r, 0.0)
        if q <= 0.0:
            print(f"    [WARNING] run {r} has zero or missing charge.")
        Q_t += q

    # 7) Mask events so only events from selected_runs remain
    mask_runs = np.isin(run_t, list(selected_runs))
    mx2_use   = mx2_t[mask_runs]

    if Q_t > 0.0 and mx2_use.size > 0:
        counts, _ = np.histogram(mx2_use, bins=bins)
        hist_t    = counts.astype(float) / Q_t

    return (label, hist_t)


# -----------------------------------------------------------------------------
# PLOTTING FUNCTION
# -----------------------------------------------------------------------------

def make_four_panel_plots(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 2×2 figure:
      - Panel (0,0): NH₃ + H₂ → Mx² histogram with |t|<1
      - Panel (0,1):  C  + H₂ → Mx² histogram with |t|<1
      - Panel (1,0): Differences [NH₃ – C] → Mx² with |t|<1
      - Panel (1,1): [(NH₃ – scaled C)/NH₃] vs x_B for 0.75<Mx²<1.05 and |t|<1

    Uses ProcessPoolExecutor to parallelize Mx²‐histogramming of each file.
    """
    # 1) Parallel histogram Mx² with |t|<1
    tasks = []
    for fp, lbl in nh3_files + c_files + h2_files:
        tasks.append((fp, lbl, run_charges, QUICK_RUN))

    results = {}  # label → hist_t
    with ProcessPoolExecutor() as executor:
        future_to_label = {executor.submit(process_file, args): args[1] for args in tasks}
        for future in as_completed(future_to_label):
            lbl = future_to_label[future]
            try:
                label, hist_t = future.result()
                results[label] = hist_t
            except Exception as e:
                print(f"[ERROR] processing {lbl}: {e}")
                results[lbl] = np.zeros(100, dtype=float)

    # 2) Build 2×2 canvas
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=False)
    bins       = np.linspace(0.0, 1.5, 101)
    bin_centers= 0.5 * (bins[:-1] + bins[1:])
    period_names = ["Su22", "Fa22", "Sp23"]
    colors = {"Su22": "C0", "Fa22": "C1", "Sp23": "C2", "Sp19": "k"}

    # ────────────────────────────────────────────────────────────────────────────
    # Panel (0,0): NH₃ + H₂ (|t|<1)
    # ────────────────────────────────────────────────────────────────────────────
    top0_vals = []
    for fp, lbl in nh3_files:
        hist = results.get(lbl, np.zeros(100))
        color = colors[lbl.split("-")[0]]
        axes[0, 0].step(bins[:-1], hist, where="post", color=color, label=lbl)
        top0_vals.append(hist.max())
    # H₂
    h2_lbl = h2_files[0][1]
    hist_h2 = results.get(h2_lbl, np.zeros(100))
    axes[0, 0].step(bins[:-1], hist_h2, where="post", color=colors["Sp19"], linestyle="-", label=h2_lbl)
    top0_vals.append(hist_h2.max())

    axes[0, 0].set(
        title="NH₃ + H₂: Normalized Mx² (|t|<1)",
        xlabel=r"$M_x^2$ (GeV^2)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    # dynamic y‐limit
    y0 = 1.1 * max(top0_vals) if top0_vals else 0.1
    axes[0, 0].set(ylim=(0.0, y0))
    axes[0, 0].legend(loc="upper right", fontsize="small")

    # ────────────────────────────────────────────────────────────────────────────
    # Panel (0,1): C + H₂ (|t|<1)
    # ────────────────────────────────────────────────────────────────────────────
    top1_vals = []
    for fp, lbl in c_files:
        hist = results.get(lbl, np.zeros(100))
        color = colors[lbl.split("-")[0]]
        axes[0, 1].step(bins[:-1], hist, where="post", color=color, label=lbl)
        top1_vals.append(hist.max())
    axes[0, 1].step(bins[:-1], hist_h2, where="post", color=colors["Sp19"], linestyle="-", label=h2_lbl)
    top1_vals.append(hist_h2.max())

    axes[0, 1].set(
        title="C + H₂: Normalized Mx² (|t|<1)",
        xlabel=r"$M_x^2$ (GeV^2)",
        xlim=(0.0, 1.5)
    )
    y1 = 1.1 * max(top1_vals) if top1_vals else 0.1
    axes[0, 1].set(ylim=(0.0, y1))
    axes[0, 1].legend(loc="upper right", fontsize="small")

    # ────────────────────────────────────────────────────────────────────────────
    # Panel (1,0): Differences [NH₃ – C] (|t|<1)
    # ────────────────────────────────────────────────────────────────────────────
    bot0_vals = []
    for per in period_names:
        nh3_lbl = f"{per}-NH3"
        c_lbl   = f"{per}-C"
        H_n = results.get(nh3_lbl, np.zeros(100))
        H_c = results.get(c_lbl,   np.zeros(100))
        # rescale C so that ∫_{0→0.5} matches
        mask_0p5 = bin_centers < 0.5
        sum_n = H_n[mask_0p5].sum()
        sum_c = H_c[mask_0p5].sum()
        scale = (sum_n / sum_c) if sum_c > 0 else 0.0
        diff = H_n - scale * H_c
        color = colors[per]
        style = "-" if per=="Su22" else ("--" if per=="Fa22" else "-.")
        axes[1, 0].step(bins[:-1], diff, where="post", color=color, linestyle=style, label=f"{per} NH₃–C")
        bot0_vals.append(diff.max())
    # Also plot H₂ (uncut)
    axes[1, 0].step(bins[:-1], hist_h2, where="post", color=colors["Sp19"], linestyle="-", label=h2_lbl)
    bot0_vals.append(hist_h2.max())

    axes[1, 0].set(
        title="(NH₃ – C) (|t|<1), + H₂",
        xlabel=r"$M_x^2$ (GeV^2)",
        xlim=(0.0, 1.5)
    )
    y2 = 1.1 * max(bot0_vals) if bot0_vals else 0.1
    axes[1, 0].set(ylim=(0.0, y2))
    axes[1, 0].legend(loc="upper right", fontsize="small")

    # ────────────────────────────────────────────────────────────────────────────
    # Panel (1,1): Ratio (H_n – scaled H_c)/H_n vs x_B
    # ────────────────────────────────────────────────────────────────────────────
    # We need to re-open each tree and collect x_B for events with:
    #    run ≤ MAX_RUNNUM, 0.75 < Mx² < 1.05, and |t| < 1.
    x_bins = np.linspace(0.0, 0.6, 26)  # 25 bins from 0→0.6
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
    period_hist = {}   # store H_n_norm and H_c_norm for each period

    for per in period_names:
        # Find filepaths
        nh3_fp = next(fp for fp, lbl in nh3_files if lbl.startswith(per))
        c_fp   = next(fp for fp, lbl in c_files   if lbl.startswith(per))

        # NH3: collect x_B
        tree_n = uproot.open(nh3_fp)["PhysicsEvents"]
        run_n  = tree_n["runnum"].array(library="np").astype(int)
        mx2_n  = tree_n["Mx2"].array(library="np")
        e_p_n  = tree_n["e_p"].array(library="np")
        e_th_n = tree_n["e_theta"].array(library="np")
        e_ph_n = tree_n["e_phi"].array(library="np")
        p_p_n  = tree_n["p_p"].array(library="np")
        p_th_n = tree_n["p_theta"].array(library="np")
        p_ph_n = tree_n["p_phi"].array(library="np")
        x_n    = tree_n["x"].array(library="np")

        mask_base_n = (run_n <= MAX_RUNNUM) & (mx2_n > 0.75) & (mx2_n < 1.05)
        run_base_n  = run_n[mask_base_n]
        mx2_base_n  = mx2_n[mask_base_n]
        e_p_base_n  = e_p_n[mask_base_n]
        e_th_base_n = e_th_n[mask_base_n]
        e_ph_base_n = e_ph_n[mask_base_n]
        p_p_base_n  = p_p_n[mask_base_n]
        p_th_base_n = p_th_n[mask_base_n]
        p_ph_base_n = p_ph_n[mask_base_n]
        x_base_n    = x_n[mask_base_n]

        if run_base_n.size > 0:
            t_n = compute_t_array(
                run_base_n,
                e_p_base_n, e_th_base_n, e_ph_base_n,
                p_p_base_n, p_th_base_n, p_ph_base_n
            )
            mask_t_n = np.abs(t_n) < 1.0
        else:
            mask_t_n = np.array([], dtype=bool)

        run_t_n = run_base_n[mask_t_n]
        x_t_n   = x_base_n[mask_t_n]

        # Only keep events from selected runs
        if QUICK_RUN:
            seen = []
            for r in run_t_n:
                if r not in seen:
                    seen.append(r)
                if len(seen) >= 5:
                    break
            sel_runs_n = set(seen)
        else:
            sel_runs_n = set(np.unique(run_t_n))

        mask_runs_n = np.isin(run_t_n, list(sel_runs_n))
        x_final_n   = x_t_n[mask_runs_n]

        # Compute Q_NH3
        Q_n = sum(run_charges.get(r, 0.0) for r in sel_runs_n)

        # Build normalized histogram in x_B
        Hn_counts, _ = np.histogram(x_final_n, bins=x_bins)
        if Q_n > 0:
            Hn_norm = Hn_counts.astype(float) / Q_n
        else:
            Hn_norm = np.zeros(len(x_centers), dtype=float)

        # C: same procedure
        tree_c = uproot.open(c_fp)["PhysicsEvents"]
        run_c  = tree_c["runnum"].array(library="np").astype(int)
        mx2_c  = tree_c["Mx2"].array(library="np")
        e_p_c  = tree_c["e_p"].array(library="np")
        e_th_c = tree_c["e_theta"].array(library="np")
        e_ph_c = tree_c["e_phi"].array(library="np")
        p_p_c  = tree_c["p_p"].array(library="np")
        p_th_c = tree_c["p_theta"].array(library="np")
        p_ph_c = tree_c["p_phi"].array(library="np")
        x_c    = tree_c["x"].array(library="np")

        mask_base_c = (run_c <= MAX_RUNNUM) & (mx2_c > 0.75) & (mx2_c < 1.05)
        run_base_c  = run_c[mask_base_c]
        mx2_base_c  = mx2_c[mask_base_c]
        e_p_base_c  = e_p_c[mask_base_c]
        e_th_base_c = e_th_c[mask_base_c]
        e_ph_base_c = e_ph_c[mask_base_c]
        p_p_base_c  = p_p_c[mask_base_c]
        p_th_base_c = p_th_c[mask_base_c]
        p_ph_base_c = p_ph_c[mask_base_c]
        x_base_c    = x_c[mask_base_c]

        if run_base_c.size > 0:
            t_c = compute_t_array(
                run_base_c,
                e_p_base_c, e_th_base_c, e_ph_base_c,
                p_p_base_c, p_th_base_c, p_ph_base_c
            )
            mask_t_c = np.abs(t_c) < 1.0
        else:
            mask_t_c = np.array([], dtype=bool)

        run_t_c = run_base_c[mask_t_c]
        x_t_c   = x_base_c[mask_t_c]

        if QUICK_RUN:
            seen = []
            for r in run_t_c:
                if r not in seen:
                    seen.append(r)
                if len(seen) >= 5:
                    break
            sel_runs_c = set(seen)
        else:
            sel_runs_c = set(np.unique(run_t_c))

        mask_runs_c = np.isin(run_t_c, list(sel_runs_c))
        x_final_c   = x_t_c[mask_runs_c]

        Q_c = sum(run_charges.get(r, 0.0) for r in sel_runs_c)

        Hc_counts, _ = np.histogram(x_final_c, bins=x_bins)
        if Q_c > 0:
            Hc_norm = Hc_counts.astype(float) / Q_c
        else:
            Hc_norm = np.zeros(len(x_centers), dtype=float)

        period_hist[per] = (Hn_norm, Hc_norm)

    # Now compute and plot (Hn_norm – scaled Hc_norm)/Hn_norm
    for per in period_names:
        Hn_norm, Hc_norm = period_hist[per]
        mask_x_lt_0p5 = x_centers < 0.5
        sum_n = Hn_norm[mask_x_lt_0p5].sum()
        sum_c = Hc_norm[mask_x_lt_0p5].sum()
        scale = (sum_n / sum_c) if sum_c > 0 else 0.0

        ratio = np.zeros_like(x_centers)
        for i in range(len(x_centers)):
            n_i = Hn_norm[i]
            c_i = Hc_norm[i]
            if n_i > 0:
                ratio[i] = (n_i - scale * c_i) / n_i
            else:
                ratio[i] = 0.0

        color = colors[per]
        axes[1, 1].step(x_centers, ratio, where="mid", color=color, linestyle=("-" if per=="Su22" else "--" if per=="Fa22" else "-."), label=f"{per}→ ratio")

        # Overlay dilution‐factor points, if available
        if per == "Su22":
            axes[1, 1].errorbar(x_Su22, dil_Su22, yerr=dil_err_Su22, fmt="o", color=color, label="Su22 dil.")
        if per == "Fa22":
            axes[1, 1].errorbar(x_Fa22, dil_Fa22, yerr=dil_err_Fa22, fmt="s", color=color, label="Fa22 dil.")
        if per == "Sp23" and x_Sp23.size > 0:
            axes[1, 1].errorbar(x_Sp23, dil_Sp23, yerr=dil_err_Sp23, fmt="^", color=color, label="Sp23 dil.")

    axes[1, 1].set(
        title=r"$(NH₃ - C)/NH₃$ vs $x_B$  (0.75<$M_x^2$<1.05, |t|<1)",
        xlabel=r"$x_B$",
        xlim=(0.0, 0.6),
        ylim=(0.0, 1.0)
    )
    axes[1, 1].legend(loc="upper right", fontsize="small")

    # 3) Save figure
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    run_charges = parse_run_charges(CSV_PATH)
    make_four_panel_plots(NH3_FILES, C_FILES, H2_FILES, run_charges, FOUR_PANEL_OUTPUT)


if __name__ == "__main__":
    main()