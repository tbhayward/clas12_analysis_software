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
THREE_PANEL_OUTPUT = os.path.join(OUTPUT_DIR, "three_panel_Mx2_comparison.pdf")

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
# UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

def parse_run_charges(csv_path):
    """
    Read clas12_run_info.csv and return a dict { runnum: charge }, skipping run > MAX_RUNNUM.
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
    Compute Mandelstam t event-by-event for ep -> e' n π⁺, assuming exclusivity.
    t = (q - p_π)², where q = p_beam - p_e'.

    - run_arr:  array of run numbers
    - e_p_arr:  array of scattered-electron momentum magnitudes (GeV)
    - e_th_arr: array of electron polar angles (radians)
    - e_ph_arr: array of electron azimuthal angles (radians)
    - p_p_arr:  array of pion momentum magnitudes (GeV)
    - p_th_arr: array of pion polar angles (radians)
    - p_ph_arr: array of pion azimuthal angles (radians)

    Returns t_array (in GeV²) of same length.
    """
    # 1) Determine Eb per event
    Eb_arr = get_beam_energy(run_arr)

    # 2) Scattered electron e' 4-vector
    E_e   = np.sqrt(e_p_arr**2 + m_e**2)
    sin_e = np.sin(e_th_arr)
    cos_e = np.cos(e_th_arr)
    ex =  e_p_arr * sin_e * np.cos(e_ph_arr)
    ey =  e_p_arr * sin_e * np.sin(e_ph_arr)
    ez =  e_p_arr * cos_e

    # 3) Pion 4-vector
    E_pi = np.sqrt(p_p_arr**2 + m_pi**2)
    sin_p = np.sin(p_th_arr)
    cos_p = np.cos(p_th_arr)
    px =  p_p_arr * sin_p * np.cos(p_ph_arr)
    py =  p_p_arr * sin_p * np.sin(p_ph_arr)
    pz =  p_p_arr * cos_p

    # 4) Virtual photon q 4-vector: q = p_beam - p_e'
    E_q = Eb_arr - E_e
    qx  = -ex
    qy  = -ey
    qz  = Eb_arr - ez

    # 5) Compute Δ = q - p_π
    dE = E_q - E_pi
    dx = qx - px
    dy = qy - py
    dz = qz - pz

    # 6) t = (ΔE)² - (Δp)²
    t_array = dE**2 - (dx**2 + dy**2 + dz**2)
    return t_array


def process_file(args):
    """
    Worker for parallel histogramming:
    Returns (label, hist_no_t, hist_t) where:
      - hist_no_t: normalized Mx² histogram with only (runnum ≤ MAX_RUNNUM, Mx² ≤ 2.0).
      - hist_t:    normalized Mx² histogram with additional |t| < 1.0 cut.
    Both histograms have 100 bins from 0 → 1.5.

    args = (filepath, label, run_charges, quick_flag)
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

    # Prepare output histograms (100 bins)
    bins = np.linspace(0.0, 1.5, 101)

    # Initialize with zeros
    hist_no_t = np.zeros(len(bins)-1, dtype=float)
    hist_t    = np.zeros(len(bins)-1, dtype=float)

    # If no events pass the base filter, return zeros
    if run_base.size == 0:
        return (label, hist_no_t, hist_t)

    # --- Compute no-t-cut histogram ---

    # 4a) Determine which runs to include (no-t): run_final_no_t = run_base
    run_final_no_t = run_base.copy()
    mx2_final_no_t = mx2_base.copy()

    # 5a) Select unique runs (first 5 if quick)
    if quick:
        seen = []
        for r in run_final_no_t:
            if r not in seen:
                seen.append(r)
            if len(seen) >= 5:
                break
        selected_runs_no_t = set(seen)
    else:
        selected_runs_no_t = set(np.unique(run_final_no_t))

    # 6a) Sum total charge Q_no_t
    Q_no_t = 0.0
    for r in selected_runs_no_t:
        q = run_charges.get(r, 0.0)
        if q <= 0.0:
            print(f"    [WARNING] (no-t) run {r} has zero or missing charge.")
        Q_no_t += q

    # 7a) Mask events to only those runs
    mask_runs_no_t = np.isin(run_final_no_t, list(selected_runs_no_t))
    mx2_use_no_t   = mx2_final_no_t[mask_runs_no_t]

    if Q_no_t > 0.0 and mx2_use_no_t.size > 0:
        counts_no_t, _ = np.histogram(mx2_use_no_t, bins=bins)
        hist_no_t = counts_no_t.astype(float) / Q_no_t

    # --- Compute t-cut histogram ---

    # 4b) Compute t-array for base-filtered events
    t_vals = compute_t_array(
        run_base,
        e_p_base, e_th_base, e_ph_base,
        p_p_base, p_th_base, p_ph_base
    )
    mask_t = np.abs(t_vals) < 1.0

    run_t     = run_base[mask_t]
    mx2_t     = mx2_base[mask_t]

    if run_t.size == 0:
        return (label, hist_no_t, hist_t)

    # 5b) Select unique runs for t-cut
    if quick:
        seen = []
        for r in run_t:
            if r not in seen:
                seen.append(r)
            if len(seen) >= 5:
                break
        selected_runs_t = set(seen)
    else:
        selected_runs_t = set(np.unique(run_t))

    # 6b) Sum total charge Q_t
    Q_t = 0.0
    for r in selected_runs_t:
        q = run_charges.get(r, 0.0)
        if q <= 0.0:
            print(f"    [WARNING] (t) run {r} has zero or missing charge.")
        Q_t += q

    # 7b) Mask events to only those runs
    mask_runs_t = np.isin(run_t, list(selected_runs_t))
    mx2_use_t   = mx2_t[mask_runs_t]

    if Q_t > 0.0 and mx2_use_t.size > 0:
        counts_t, _ = np.histogram(mx2_use_t, bins=bins)
        hist_t = counts_t.astype(float) / Q_t

    return (label, hist_no_t, hist_t)


# -----------------------------------------------------------------------------
# PLOTTING FUNCTION
# -----------------------------------------------------------------------------

def make_normalized_Mx2_plots(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 2×3 figure (two rows, three columns):
      - Top row: no |t| < 1 cut (only runnum ≤ MAX_RUNNUM, Mx² ≤ 2.0)
        • Panel (0,0): NH₃ + H₂
        • Panel (0,1):  C   + H₂
        • Panel (0,2): Differences [NH₃ – C] + H₂
      - Bottom row: with |t| < 1 cut (in addition to above)
        • Panel (1,0): NH₃ + H₂ (|t|<1)
        • Panel (1,1):  C   + H₂ (|t|<1)
        • Panel (1,2): Differences [NH₃ – C] + H₂ (|t|<1)

    Uses ProcessPoolExecutor to parallelize histogramming of each file.
    """
    # 1) Prepare tasks for parallel execution
    tasks = []
    for fp, lbl in nh3_files + c_files + h2_files:
        tasks.append((fp, lbl, run_charges, QUICK_RUN))

    # 2) Submit tasks
    results = {}  # label → (hist_no_t, hist_t)
    with ProcessPoolExecutor() as executor:
        future_to_label = {executor.submit(process_file, args): args[1] for args in tasks}
        for future in as_completed(future_to_label):
            lbl = future_to_label[future]
            try:
                label, hist_no_t, hist_t = future.result()
                results[label] = (hist_no_t, hist_t)
            except Exception as e:
                print(f"[ERROR] processing {lbl}: {e}")
                # In case of error, store zero arrays
                results[lbl] = (np.zeros(100, dtype=float), np.zeros(100, dtype=float))

    # 3) Prepare bins and bin centers
    bins = np.linspace(0.0, 1.5, 101)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    # 4) Start plotting: 2 rows × 3 columns
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharey=True)

    # Prepare lists of all histograms for each row
    # Top row: collect NH₃ no-t, C no-t, diff no-t, H₂ no-t
    top_row_vals = []
    # Bottom row: collect NH₃ with-t, C with-t, diff with-t, H₂ with-t
    bottom_row_vals = []

    # ===========================
    # Top row: no |t| < 1 cut
    # ===========================

    # Panel (0,0): NH₃ + H₂   (no-t)
    for _, lbl in nh3_files:
        hist_no_t, _ = results.get(lbl, (np.zeros(100), np.zeros(100)))
        axes[0, 0].step(bins[:-1], hist_no_t, where='post', label=lbl)
        top_row_vals.append(hist_no_t.max())
    h2_label = h2_files[0][1]
    hist_h2_no_t, _ = results.get(h2_label, (np.zeros(100), np.zeros(100)))
    axes[0, 0].step(
        bins[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    top_row_vals.append(hist_h2_no_t.max())
    axes[0, 0].set(
        title="NH₃ + H₂: Normalized Mx² (no t-cut)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    axes[0, 0].legend(loc='upper right', fontsize='small')

    # Panel (0,1): C + H₂   (no-t)
    for _, lbl in c_files:
        hist_no_t, _ = results.get(lbl, (np.zeros(100), np.zeros(100)))
        axes[0, 1].step(bins[:-1], hist_no_t, where='post', label=lbl)
        top_row_vals.append(hist_no_t.max())
    axes[0, 1].step(
        bins[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    top_row_vals.append(hist_h2_no_t.max())
    axes[0, 1].set(
        title="C + H₂: Normalized Mx² (no t-cut)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    axes[0, 1].legend(loc='upper right', fontsize='small')

    # Panel (0,2): Differences [NH₃ – C] + H₂   (no-t)
    period_names = ["Su22", "Fa22", "Sp23"]
    nh3_hist_no_t = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[0] for _, lbl in nh3_files}
    c_hist_no_t   = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[0] for _, lbl in c_files}

    diff_no_t = {}
    for period_label in period_names:
        nh3_lbl = f"{period_label}-NH3"
        c_lbl   = f"{period_label}-C"
        nh3_vals = nh3_hist_no_t.get(nh3_lbl, np.zeros(100))
        c_vals   = c_hist_no_t.get(c_lbl,   np.zeros(100))

        mask_0_5 = (bin_centers >= 0.0) & (bin_centers < 0.5)
        nh3_sum = nh3_vals[mask_0_5].sum()
        c_sum   = c_vals[mask_0_5].sum()
        scale = (nh3_sum / c_sum) if c_sum > 0.0 else 0.0
        print(f"[DEBUG] (no-t) {period_label}: NH₃ sum(0–0.5)={nh3_sum:.3f}, C sum(0–0.5)={c_sum:.3f}, scale={scale:.3f}")
        diff_no_t[period_label] = nh3_vals - (c_vals * scale)

    for idx, period_label in enumerate(period_names):
        diff = diff_no_t.get(period_label, np.zeros(100))
        style = ['-', '--', '-.', ':'][idx % 4]
        axes[0, 2].step(
            bins[:-1], diff, where='post',
            linestyle=style, label=f"{period_label} NH₃–C"
        )
        top_row_vals.append(diff.max())

    axes[0, 2].step(
        bins[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    top_row_vals.append(hist_h2_no_t.max())

    # Determine top row y-limit = 1.1 × maximum of top_row_vals
    y_top = 1.1 * max(top_row_vals) if top_row_vals else 0.1
    for col in range(3):
        axes[0, col].set(ylim=(0.0, y_top))

    axes[0, 2].set(
        title="(NH₃ – C) + H₂ (no t-cut)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    axes[0, 2].legend(loc='upper right', fontsize='small')

    # ===========================
    # Bottom row: with |t| < 1 cut
    # ===========================

    # Reset lists for bottom row maxima
    bottom_row_vals = []

    # Panel (1,0): NH₃ + H₂   (with-t)
    for _, lbl in nh3_files:
        _, hist_t = results.get(lbl, (np.zeros(100), np.zeros(100)))
        axes[1, 0].step(bins[:-1], hist_t, where='post', label=lbl)
        bottom_row_vals.append(hist_t.max())
    _, hist_h2_t = results.get(h2_label, (np.zeros(100), np.zeros(100)))
    axes[1, 0].step(
        bins[:-1], hist_h2_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    bottom_row_vals.append(hist_h2_t.max())
    axes[1, 0].set(
        title="NH₃ + H₂: Normalized Mx² (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    axes[1, 0].legend(loc='upper right', fontsize='small')

    # Panel (1,1): C + H₂   (with-t)
    for _, lbl in c_files:
        _, hist_t = results.get(lbl, (np.zeros(100), np.zeros(100)))
        axes[1, 1].step(bins[:-1], hist_t, where='post', label=lbl)
        bottom_row_vals.append(hist_t.max())
    axes[1, 1].step(
        bins[:-1], hist_h2_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    bottom_row_vals.append(hist_h2_t.max())
    axes[1, 1].set(
        title="C + H₂: Normalized Mx² (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    axes[1, 1].legend(loc='upper right', fontsize='small')

    # Panel (1,2): Differences [NH₃ – C] + H₂   (with-t)
    nh3_hist_t = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[1] for _, lbl in nh3_files}
    c_hist_t   = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[1] for _, lbl in c_files}

    diff_t = {}
    for period_label in period_names:
        nh3_lbl = f"{period_label}-NH3"
        c_lbl   = f"{period_label}-C"
        nh3_vals = nh3_hist_t.get(nh3_lbl, np.zeros(100))
        c_vals   = c_hist_t.get(c_lbl,   np.zeros(100))

        mask_0_5 = (bin_centers >= 0.0) & (bin_centers < 0.5)
        nh3_sum = nh3_vals[mask_0_5].sum()
        c_sum   = c_vals[mask_0_5].sum()
        scale = (nh3_sum / c_sum) if c_sum > 0.0 else 0.0
        print(f"[DEBUG] (with-t) {period_label}: NH₃ sum(0–0.5)={nh3_sum:.3f}, C sum(0–0.5)={c_sum:.3f}, scale={scale:.3f}")
        diff_t[period_label] = nh3_vals - (c_vals * scale)

    for idx, period_label in enumerate(period_names):
        diff = diff_t.get(period_label, np.zeros(100))
        style = ['-', '--', '-.', ':'][idx % 4]
        axes[1, 2].step(
            bins[:-1], diff, where='post',
            linestyle=style, label=f"{period_label} NH₃–C"
        )
        bottom_row_vals.append(diff.max())
    axes[1, 2].step(
        bins[:-1], hist_h2_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    bottom_row_vals.append(hist_h2_t.max())

    # Determine bottom row y-limit = 1.1 × maximum of bottom_row_vals
    y_bottom = 1.1 * max(bottom_row_vals) if bottom_row_vals else 0.1
    for col in range(3):
        axes[1, col].set(ylim=(0.0, y_bottom))

    axes[1, 2].set(
        title="(NH₃ – C) + H₂ (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    axes[1, 2].legend(loc='upper right', fontsize='small')

    # 5) Save the 2×3 figure
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    # Load run charges
    run_charges = parse_run_charges(CSV_PATH)

    # Create the 2×3 comparison plot (parallelized)
    make_normalized_Mx2_plots(NH3_FILES, C_FILES, H2_FILES, run_charges, THREE_PANEL_OUTPUT)

if __name__ == "__main__":
    main()