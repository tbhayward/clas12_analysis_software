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
QUICK_RUN = False

# Maximum run number to include
MAX_RUNNUM = 17768

# Path to the CSV of run charges
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output directory and filename
OUTPUT_DIR            = "output/enpi+"
THREE_PANEL_OUTPUT    = os.path.join(OUTPUT_DIR, "three_panel_Mx2_comparison.pdf")
PERIOD_COMPARISON_OUT = os.path.join(OUTPUT_DIR, "period_NH3_C_Mx2_comparison.pdf")
XB_RATIO_OUTPUT       = os.path.join(OUTPUT_DIR, "xB_ratio_comparison.pdf")

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
# ASYMMETRY FIT DATA & DILUTION FACTORS
# -----------------------------------------------------------------------------

# RGC Su22 data (new iteration)
enpichi2FitsALUsinphi_Su22 = [
    [0.094266731, 1.184739999, 0.126228827],
    [0.168622194, 0.082358508, 0.012946927],
    [0.254887003, 0.115297657, 0.008854289],
    [0.348247742, 0.144713802, 0.009273995],
    [0.441277333, 0.117431442, 0.012353808],
    [0.535009222, 0.094723130, 0.021626204]
]

# RGC Fa22 data (updated results)
enpichi2FitsALUsinphi_Fa22 = [
    [0.094825648, 1.221850467, 0.095501230],
    [0.168249360, 0.081525475, 0.009989312],
    [0.255099359, 0.118300134, 0.007065332],
    [0.348334763, 0.140055462, 0.007378878],
    [0.441201661, 0.140928523, 0.009773392],
    [0.534977263, 0.110737425, 0.016852914]
]

# RGC Sp23 data (newly provided)
enpichi2FitsALUsinphi_Sp23 = [
    [0.091314751, 0.060825979, 0.081001615],
    [0.166606555, 0.070263513, 0.012697034],
    [0.251219557, 0.107175756, 0.009614248],
    [0.346388218, 0.141913308, 0.010857952],
    [0.441235552, 0.118073067, 0.014919450],
    [0.535187978, 0.115828125, 0.026600146]
]

# Dilution factor data for Su22
x_Su22 = np.array([row[0] for row in enpichi2FitsALUsinphi_Su22])
dil_Su22 = np.array([0.711119, 0.376873, 0.390790, 0.401177, 0.410563, 0.416077])
dil_err_Su22 = np.array([0.173083, 0.0175818, 0.00942491, 0.00801711, 0.0108668, 0.0225813])

# Dilution factor data for Fa22
x_Fa22 = np.array([row[0] for row in enpichi2FitsALUsinphi_Fa22])
dil_Fa22 = np.array([0.374532, 0.406101, 0.386122, 0.397153, 0.415526, 0.428239])
dil_err_Fa22 = np.array([0.117283, 0.00608866, 0.00354573, 0.00301186, 0.00400862, 0.00824607])

# For Sp23, leave blank arrays (will overlay when available)
x_Sp23 = np.array([])
dil_Sp23 = np.array([])
dil_err_Sp23 = np.array([])

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

    # 3) Virtual photon q 4-vector: q = p_beam – p_e'
    E_q = Eb_arr - E_e
    qx  = -ex
    qy  = -ey
    qz  = Eb_arr - ez

    # 4) Δ = q – p_π
    dE = E_q - E_pi
    dx = qx - px
    dy = qy - py
    dz = qz - pz

    # 5) t = (ΔE)² – (Δp)²
    t_array = dE**2 - (dx**2 + dy**2 + dz**2)
    return t_array


# -----------------------------------------------------------------------------
# PROCESSING FUNCTIONS
# -----------------------------------------------------------------------------

def process_file(args):
    """
    Worker for parallel histogramming:
    Returns (label, hist_no_t, hist_t) where:
      - hist_no_t: normalized Mx² histogram with just (runnum ≤ MAX_RUNNUM, Mx² ≤ 2.0).
      - hist_t:    normalized Mx² histogram with additional |t| < 1.0 cut.
    Both histograms use 100 bins from 0 → 1.5.

    args = (filepath, label, run_charges, quick_flag)
    """
    filepath, label, run_charges, quick = args

    if quick:
        return process_file_quick_stop(filepath, label, run_charges)

    # -----------------------------------------------------------------------------
    # Full-array version (quick=False)
    # -----------------------------------------------------------------------------

    # 1) Open TTree
    tree = uproot.open(filepath)["PhysicsEvents"]

    # 2) Extract entire arrays from the TTree
    run_arr  = tree["runnum"].array(library="np").astype(int)
    mx2_arr  = tree["Mx2"].array(library="np")
    e_p_arr  = tree["e_p"].array(library="np")
    e_th_arr = tree["e_theta"].array(library="np")
    e_ph_arr = tree["e_phi"].array(library="np")
    p_p_arr  = tree["p_p"].array(library="np")
    p_th_arr = tree["p_theta"].array(library="np")
    p_ph_arr = tree["p_phi"].array(library="np")

    # 3) Base filter: keep only events with run ≤ MAX_RUNNUM and Mx² ≤ 2.0
    mask_base = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
    run_base  = run_arr[mask_base]
    mx2_base  = mx2_arr[mask_base]
    e_p_base  = e_p_arr[mask_base]
    e_th_base = e_th_arr[mask_base]
    e_ph_base = e_ph_arr[mask_base]
    p_p_base  = p_p_arr[mask_base]
    p_th_base = p_th_arr[mask_base]
    p_ph_base = p_ph_arr[mask_base]

    # Print how many events remain after the base filter
    num_base_events = run_base.size
    print(f"[DEBUG][{label}] after base filter: {num_base_events} events")

    # Prepare output histograms (100 bins: 0.0 → 1.5)
    bins = np.linspace(0.0, 1.5, 101)
    hist_no_t = np.zeros(len(bins) - 1, dtype=float)
    hist_t    = np.zeros(len(bins) - 1, dtype=float)

    # If no events survive the base filter, return zeroed histograms
    if num_base_events == 0:
        print(f"[DEBUG][{label}] no events remain after base filter → all histograms zero")
        return (label, hist_no_t, hist_t)

    # --- Compute “no‐t‐cut” histogram (only run/Mx² filter) ---
    run_final_no_t = run_base.copy()
    mx2_final_no_t = mx2_base.copy()

    # 5a) Select which runs contribute: all unique runs (quick=False)
    selected_runs_no_t = set(np.unique(run_final_no_t))

    # 6a) Sum total beam charge for those runs (Q_no_t) – still needed for normalization
    Q_no_t = 0.0
    for r in selected_runs_no_t:
        q = run_charges.get(r, 0.0)
        if q <= 0.0:
            print(f"    [WARNING] (no‐t) run {r} has zero or missing charge.")
        Q_no_t += q

    # 7a) Mask events so that only those runs remain
    mask_runs_no_t = np.isin(run_final_no_t, list(selected_runs_no_t))
    mx2_use_no_t   = mx2_final_no_t[mask_runs_no_t]

    if Q_no_t > 0.0 and mx2_use_no_t.size > 0:
        counts_no_t, _ = np.histogram(mx2_use_no_t, bins=bins)
        hist_no_t      = counts_no_t.astype(float) / Q_no_t
    else:
        hist_no_t = np.zeros(len(bins) - 1, dtype=float)

    # --- Compute “t‐cut” histogram (also require |t|<1) ---
    t_vals = compute_t_array(
        run_base,
        e_p_base, e_th_base, e_ph_base,
        p_p_base, p_th_base, p_ph_base
    )

    mask_t = np.abs(t_vals) < 1.0
    run_t  = run_base[mask_t]
    mx2_t  = mx2_base[mask_t]

    # Count how many survive |t|<1
    num_t_events = run_t.size
    print(f"[DEBUG][{label}] after |t|<1 cut: {num_t_events} events")

    if num_t_events != 0:
        # 5b) Select which runs now contribute under |t|<1
        selected_runs_t = set(np.unique(run_t))

        # 6b) Sum total beam charge for those runs (Q_t)
        Q_t = 0.0
        for r in selected_runs_t:
            q = run_charges.get(r, 0.0)
            if q <= 0.0:
                print(f"    [WARNING] (t) run {r} has zero or missing charge.")
            Q_t += q

        # 7b) Mask events so that only those runs remain
        mask_runs_t = np.isin(run_t, list(selected_runs_t))
        mx2_use_t   = mx2_t[mask_runs_t]

        if Q_t > 0.0 and mx2_use_t.size > 0:
            counts_t, _ = np.histogram(mx2_use_t, bins=bins)
            hist_t      = counts_t.astype(float) / Q_t
        else:
            hist_t = np.zeros(len(bins) - 1, dtype=float)

    # Compute survival fraction in terms of event count
    survival_frac = 100.0 * num_t_events / num_base_events if num_base_events > 0 else 0.0
    print(f"[DEBUG][{label}] event survival = {survival_frac:.1f}% ({num_t_events}/{num_base_events})")

    return (label, hist_no_t, hist_t)


def process_file_quick_stop(filepath, label, run_charges):
    """
    Quick‐stop version of process_file():
    - Only process the FIRST 5 unique runs to build both histograms.
    - Uses uproot.iterate() to read chunks of events, stops as soon as 5 distinct runs are found.
    - Returns (label, hist_no_t, hist_t).
    """
    bins = np.linspace(0.0, 1.5, 101)
    hist_no_t = np.zeros(len(bins) - 1, dtype=float)
    hist_t    = np.zeros(len(bins) - 1, dtype=float)

    seen_runs_no_t = set()
    seen_runs_t    = set()
    mx2_acc_no_t   = []
    mx2_acc_t      = []

    total_base_events = 0

    # Iterate in chunks of 100k events, stop once each set has 5 distinct runs
    for chunk in uproot.iterate(
            filepath + ":PhysicsEvents",
            ["runnum", "Mx2", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi"],
            step_size=100_000
        ):
        run_chunk  = np.array(chunk["runnum"], dtype=int)
        mx2_chunk  = np.array(chunk["Mx2"], dtype=float)
        e_p_chunk  = np.array(chunk["e_p"], dtype=float)
        e_th_chunk = np.array(chunk["e_theta"], dtype=float)
        e_ph_chunk = np.array(chunk["e_phi"], dtype=float)
        p_p_chunk  = np.array(chunk["p_p"], dtype=float)
        p_th_chunk = np.array(chunk["p_theta"], dtype=float)
        p_ph_chunk = np.array(chunk["p_phi"], dtype=float)

        # Base filter: run ≤ MAX_RUNNUM & Mx² ≤ 2.0
        mask_base = (run_chunk <= MAX_RUNNUM) & (mx2_chunk <= 2.0)
        run_base_chunk  = run_chunk[mask_base]
        mx2_base_chunk  = mx2_chunk[mask_base]
        e_p_base_chunk  = e_p_chunk[mask_base]
        e_th_base_chunk = e_th_chunk[mask_base]
        e_ph_base_chunk = e_ph_chunk[mask_base]
        p_p_base_chunk  = p_p_chunk[mask_base]
        p_th_base_chunk = p_th_chunk[mask_base]
        p_ph_base_chunk = p_ph_chunk[mask_base]

        total_base_events += run_base_chunk.size

        # Loop event-by-event until we have 5 distinct runs in each category
        for (rnum, mx2_val, ep, eth, eph, pp, pth, pph) in zip(
                run_base_chunk, mx2_base_chunk,
                e_p_base_chunk, e_th_base_chunk, e_ph_base_chunk,
                p_p_base_chunk, p_th_base_chunk, p_ph_base_chunk
            ):

            # (1) “no-t” side: if we’ve seen <5 runs, possibly add run & Mx²
            if len(seen_runs_no_t) < 5:
                if rnum not in seen_runs_no_t:
                    seen_runs_no_t.add(rnum)
                if rnum in seen_runs_no_t:
                    mx2_acc_no_t.append(mx2_val)

            # (2) Compute t for this event only if we still need runs for the t-cut
            if len(seen_runs_t) < 5:
                t_val = compute_t_array(
                    np.array([rnum]),
                    np.array([ep]), np.array([eth]), np.array([eph]),
                    np.array([pp]), np.array([pth]), np.array([pph])
                )[0]
                if abs(t_val) < 1.0:
                    if rnum not in seen_runs_t:
                        seen_runs_t.add(rnum)
                    if rnum in seen_runs_t:
                        mx2_acc_t.append(mx2_val)

            # If both categories already have 5 runs, stop reading further events
            if len(seen_runs_no_t) >= 5 and len(seen_runs_t) >= 5:
                break

        if len(seen_runs_no_t) >= 5 and len(seen_runs_t) >= 5:
            break

    # Now we know how many base‐filtered events we saw (though we only recorded Mx² from the first 5 runs)
    num_base_events = total_base_events
    num_no_t_events = len(mx2_acc_no_t)
    num_t_events    = len(mx2_acc_t)

    # Build “no‐t” histogram
    seen_runs_list_no_t = list(seen_runs_no_t)
    Q_no_t = sum(run_charges.get(r, 0.0) for r in seen_runs_list_no_t)
    if Q_no_t > 0.0 and num_no_t_events > 0:
        counts_no_t, _ = np.histogram(np.array(mx2_acc_no_t), bins=bins)
        hist_no_t      = counts_no_t.astype(float) / Q_no_t
    else:
        hist_no_t = np.zeros(len(bins) - 1, dtype=float)

    # Build “t‐cut” histogram
    seen_runs_list_t = list(seen_runs_t)
    Q_t = sum(run_charges.get(r, 0.0) for r in seen_runs_list_t)
    if Q_t > 0.0 and num_t_events > 0:
        counts_t, _ = np.histogram(np.array(mx2_acc_t), bins=bins)
        hist_t      = counts_t.astype(float) / Q_t
    else:
        hist_t = np.zeros(len(bins) - 1, dtype=float)

    # Print out counts and survival fraction
    survival_frac = 100.0 * num_t_events / num_base_events if num_base_events > 0 else 0.0
    print(f"[DEBUG][{label}, quick] base_e = {num_base_events}  |  t_e = {num_t_events}  |  "
          f"survival = {survival_frac:.1f}% ({num_t_events}/{num_base_events})")

    return (label, hist_no_t, hist_t)


# -----------------------------------------------------------------------------
# PLOTTING FUNCTION
# -----------------------------------------------------------------------------

def make_normalized_Mx2_plots(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 2×3 figure (two rows, three columns):
      - Top row: no |t| < 1 cut           (runnum ≤ MAX_RUNNUM, Mx² ≤ 2.0)
        • Panel(0,0): NH₃ + H₂
        • Panel(0,1):  C   + H₂
        • Panel(0,2): Differences [NH₃–s·C] + H₂
      - Bottom row: with |t| < 1 cut     (in addition to above)
        • Panel(1,0): NH₃ + H₂ (|t|<1)
        • Panel(1,1):  C   + H₂ (|t|<1)
        • Panel(1,2): Differences [NH₃–s·C] + H₂ (|t|<1)

    Uses ProcessPoolExecutor to parallelize histogramming of each file,
    then also makes a separate 1×3 figure comparing NH₃ vs C for each period (no t-cut),
    and finally a single-panel xB–ratio plot with dilution factors overlay.
    """
    # 1) Prepare tasks for Mx² histograms
    tasks = []
    for fp, lbl in nh3_files + c_files + h2_files:
        tasks.append((fp, lbl, run_charges, QUICK_RUN))

    # 2) Submit tasks in parallel
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
                results[lbl] = (np.zeros(100, dtype=float), np.zeros(100, dtype=float))

    # 3) Prepare bins & bin centers for Mx²
    bins     = np.linspace(0.0, 1.5, 101)
    centers  = 0.5 * (bins[:-1] + bins[1:])

    # 4) Create the 2×3 figure (no shared y between rows)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharey=False)

    # Keep track of maxima in each row to set dynamic y‐limits
    top_row_vals    = []
    bottom_row_vals = []

    ## ────────────────────────────────────────────────────────────────────────────
    # Top row: no |t| < 1 cut
    ## ────────────────────────────────────────────────────────────────────────────

    # Panel (0,0): NH₃ + H₂ (no-t)
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
        title="NH₃ + H₂: Normalized $M_x^2$ (no t-cut)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    axes[0, 0].legend(loc='upper right', fontsize='small')

    # Panel (0,1): C + H₂ (no-t)
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
        title="C + H₂: Normalized $M_x^2$ (no t-cut)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    axes[0, 1].legend(loc='upper right', fontsize='small')

    # Panel (0,2): Differences [NH₃ – s·C] + H₂ (no-t)
    period_names    = ["Su22", "Fa22", "Sp23"]
    nh3_hist_no_t   = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[0] for _, lbl in nh3_files}
    c_hist_no_t     = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[0] for _, lbl in c_files}
    diff_no_t       = {}
    scale_dict_no   = {}

    for period_label in period_names:
        nh3_lbl  = f"{period_label}-NH3"
        c_lbl    = f"{period_label}-C"
        nh3_vals = nh3_hist_no_t.get(nh3_lbl, np.zeros(100))
        c_vals   = c_hist_no_t.get(c_lbl,   np.zeros(100))

        mask_0_5 = (centers >= 0.0) & (centers < 0.5)
        nh3_sum  = nh3_vals[mask_0_5].sum()
        c_sum    = c_vals[mask_0_5].sum()
        scale    = (nh3_sum / c_sum) if c_sum > 0.0 else 0.0
        scale_dict_no[period_label] = scale
        print(f"[DEBUG] (no-t) {period_label}: NH₃ sum(0–0.5)={nh3_sum:.3f}, C sum(0–0.5)={c_sum:.3f}, scale={scale:.3f}")
        diff_no_t[period_label] = nh3_vals - (c_vals * scale)

    for idx, period_label in enumerate(period_names):
        diff  = diff_no_t.get(period_label, np.zeros(100))
        scale = scale_dict_no[period_label]
        style = ['-', '--', '-.', ':'][idx % 4]
        axes[0, 2].step(
            bins[:-1], diff, where='post',
            linestyle=style,
            label=f"{period_label} NH₃ – {scale:.3f}·C"
        )
        top_row_vals.append(diff.max())

    axes[0, 2].step(
        bins[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    top_row_vals.append(hist_h2_no_t.max())

    # Determine top‐row y‐limit = 1.1 × maximum of top_row_vals
    y_top = 1.1 * max(top_row_vals) if top_row_vals else 0.1
    for col in range(3):
        axes[0, col].set(ylim=(0.0, y_top))

    axes[0, 2].set(
        title="(NH₃ – s·C) + H₂ (no t-cut)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    axes[0, 2].legend(loc='upper right', fontsize='small')

    ## ────────────────────────────────────────────────────────────────────────────
    # Bottom row: with |t| < 1 cut
    ## ────────────────────────────────────────────────────────────────────────────

    bottom_row_vals = []

    # Panel (1,0): NH₃ + H₂ (|t|<1)
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
        title="NH₃ + H₂: Normalized $M_x^2$ (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    axes[1, 0].legend(loc='upper right', fontsize='small')

    # Panel (1,1): C + H₂ (|t|<1)
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
        title="C + H₂: Normalized $M_x^2$ (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    axes[1, 1].legend(loc='upper right', fontsize='small')

    # Panel (1,2): Differences [NH₃ – s·C] + H₂ (|t|<1)
    nh3_hist_t     = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[1] for _, lbl in nh3_files}
    c_hist_t       = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[1] for _, lbl in c_files}
    diff_t         = {}
    scale_dict_t   = {}

    for period_label in period_names:
        nh3_lbl  = f"{period_label}-NH3"
        c_lbl    = f"{period_label}-C"
        nh3_vals = nh3_hist_t.get(nh3_lbl, np.zeros(100))
        c_vals   = c_hist_t.get(c_lbl,   np.zeros(100))

        mask_0_5 = (centers >= 0.0) & (centers < 0.5)
        nh3_sum  = nh3_vals[mask_0_5].sum()
        c_sum    = c_vals[mask_0_5].sum()
        scale    = (nh3_sum / c_sum) if c_sum > 0.0 else 0.0
        scale_dict_t[period_label] = scale
        print(f"[DEBUG] (|t|<1) {period_label}: NH₃ sum(0–0.5)={nh3_sum:.3f}, C sum(0–0.5)={c_sum:.3f}, scale={scale:.3f}")
        diff_t[period_label] = nh3_vals - (c_vals * scale)

    for idx, period_label in enumerate(period_names):
        diff  = diff_t.get(period_label, np.zeros(100))
        scale = scale_dict_t[period_label]
        style = ['-', '--', '-.', ':'][idx % 4]
        axes[1, 2].step(
            bins[:-1], diff, where='post',
            linestyle=style,
            label=f"{period_label} NH₃ – {scale:.3f}·C"
        )
        bottom_row_vals.append(diff.max())

    axes[1, 2].step(
        bins[:-1], hist_h2_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    bottom_row_vals.append(hist_h2_t.max())

    # Determine bottom‐row y‐limit = 1.1 × maximum of bottom_row_vals
    y_bottom = 1.1 * max(bottom_row_vals) if bottom_row_vals else 0.1
    for col in range(3):
        axes[1, col].set(ylim=(0.0, y_bottom))

    axes[1, 2].set(
        title="(NH₃ – s·C) + H₂ (|t|<1)",
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
    # 1×3 “period-by-period” NH₃ vs C comparison (no t-cut)
    # -----------------------------------------------------------------------------

    bins_mx2    = bins
    centers_mx2 = centers

    fig2, axs2 = plt.subplots(1, 3, figsize=(18, 6), sharey=False)
    for idx, period_label in enumerate(period_names):
        ax = axs2[idx]
        nh3_lbl = f"{period_label}-NH3"
        c_lbl   = f"{period_label}-C"

        hist_nh3, _ = results.get(nh3_lbl, (np.zeros(100), np.zeros(100)))
        hist_c,   _ = results.get(c_lbl,   (np.zeros(100), np.zeros(100)))

        ax.step(bins_mx2[:-1], hist_nh3, where='post', label=f"{period_label}-NH₃")
        ax.step(bins_mx2[:-1], hist_c,   where='post', label=f"{period_label}-C", linestyle='--')

        ax.set(
            title=f"{period_label}: NH₃ vs C (no t-cut)",
            xlabel=r"$M_x^2$ (GeV$^2$)",
            xlim=(0.0, 1.5)
        )
        # set y-limit to encompass both histograms
        ymax = max(hist_nh3.max(), hist_c.max()) * 1.1 if (hist_nh3.size and hist_c.size) else 0.1
        ax.set_ylim(0.0, ymax)
        ax.legend(loc='upper right', fontsize='small')

    axs2[0].set_ylabel("events / nC")
    plt.tight_layout()
    plt.savefig(PERIOD_COMPARISON_OUT)
    plt.close()

    # -----------------------------------------------------------------------------
    # Single‐panel xB‐ratio plot with dilution factors overlay
    # -----------------------------------------------------------------------------

    # We'll reuse s·C scales from the |t|<1 bottom row: scale_dict_t
    # Build NH₃(x) and C(x) histograms under cuts: 0.75 < Mx² < 1.05, |t| < 1.
    x_bins = np.linspace(0.0, 1.0, 51)
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

    # Containers for normalized counts
    nh3_counts = { "Su22": np.zeros(len(x_bins)-1, dtype=float),
                   "Fa22": np.zeros(len(x_bins)-1, dtype=float),
                   "Sp23": np.zeros(len(x_bins)-1, dtype=float) }
    c_counts   = { "Su22": np.zeros(len(x_bins)-1, dtype=float),
                   "Fa22": np.zeros(len(x_bins)-1, dtype=float),
                   "Sp23": np.zeros(len(x_bins)-1, dtype=float) }
    nh3_charge = { "Su22": 0.0, "Fa22": 0.0, "Sp23": 0.0 }
    c_charge   = { "Su22": 0.0, "Fa22": 0.0, "Sp23": 0.0 }

    # Helper to fill counts for one set of files (NH3 or C)
    def fill_x_counts(files, counts_dict, charge_dict, prefix):
        for filepath, label in files:
            period = label.split('-')[0]  # e.g. "Su22"
            tree = uproot.open(filepath)["PhysicsEvents"]

            run_arr  = tree["runnum"].array(library="np").astype(int)
            mx2_arr  = tree["Mx2"].array(library="np")
            x_arr    = tree["x"].array(library="np")
            e_p_arr  = tree["e_p"].array(library="np")
            e_th_arr = tree["e_theta"].array(library="np")
            e_ph_arr = tree["e_phi"].array(library="np")
            p_p_arr  = tree["p_p"].array(library="np")
            p_th_arr = tree["p_theta"].array(library="np")
            p_ph_arr = tree["p_phi"].array(library="np")

            # 1) Cut: run ≤ MAX_RUNNUM, 0.75 < Mx² < 1.05
            mask_base_x = (run_arr <= MAX_RUNNUM) & (mx2_arr > 0.75) & (mx2_arr < 1.05)
            run_bx  = run_arr[mask_base_x]
            mx2_bx  = mx2_arr[mask_base_x]
            x_bx    = x_arr[mask_base_x]
            e_p_bx  = e_p_arr[mask_base_x]
            e_th_bx = e_th_arr[mask_base_x]
            e_ph_bx = e_ph_arr[mask_base_x]
            p_p_bx  = p_p_arr[mask_base_x]
            p_th_bx = p_th_arr[mask_base_x]
            p_ph_bx = p_ph_arr[mask_base_x]

            if run_bx.size == 0:
                continue

            # 2) Compute t-array, then cut |t|<1
            t_vals_x = compute_t_array(
                run_bx,
                e_p_bx, e_th_bx, e_ph_bx,
                p_p_bx, p_th_bx, p_ph_bx
            )
            mask_t_x = np.abs(t_vals_x) < 1.0
            run_tx = run_bx[mask_t_x]
            x_tx   = x_bx[mask_t_x]

            if run_tx.size == 0:
                continue

            # 3) Determine unique runs, sum charges
            unique_runs = np.unique(run_tx)
            for r in unique_runs:
                q = run_charges.get(r, 0.0)
                if q <= 0.0:
                    print(f"    [WARNING][{prefix},{period}] run {r} has zero or missing charge.")
                charge_dict[period] += q

            # 4) Fill x histogram only for events whose run ∈ unique_runs
            mask_select = np.isin(run_tx, unique_runs)
            x_select = x_tx[mask_select]
            if x_select.size > 0:
                counts_x, _ = np.histogram(x_select, bins=x_bins)
                counts_dict[period] += counts_x

    # Fill NH₃ counts
    fill_x_counts(NH3_FILES, nh3_counts, nh3_charge, "NH₃")
    # Fill C counts
    fill_x_counts(C_FILES, c_counts, c_charge, "C")

    # 5) Build normalized NH₃(x) and C(x)
    nh3_norm = {}
    c_norm   = {}
    for period in period_names:
        nh3_norm[period] = np.zeros_like(nh3_counts[period], dtype=float)
        c_norm[period]   = np.zeros_like(c_counts[period], dtype=float)
        if nh3_charge[period] > 0.0:
            nh3_norm[period] = nh3_counts[period].astype(float) / nh3_charge[period]
        if c_charge[period] > 0.0:
            c_norm[period] = c_counts[period].astype(float) / c_charge[period]

    # 6) Compute (NH₃(x) - s·C(x)) / NH₃(x) for each period, using scale from scale_dict_t
    ratio_x = {}
    for period in period_names:
        ratio_x[period] = np.zeros_like(nh3_norm[period], dtype=float)
        s = scale_dict_t.get(period, 0.0)
        for i in range(len(ratio_x[period])):
            nh3_val = nh3_norm[period][i]
            c_val   = c_norm[period][i]
            if nh3_val > 0.0:
                ratio_x[period][i] = (nh3_val - s * c_val) / nh3_val
            else:
                ratio_x[period][i] = 0.0

    # 7) Plot single‐panel xB‐ratio
    fig3, ax3 = plt.subplots(1, 1, figsize=(8, 6))

    # We will plot lines for each period, capture their colors
    line_handles = {}
    for period in period_names:
        line, = ax3.step(
            x_centers, ratio_x[period], where='mid',
            label=f"{period} (NH₃ – {scale_dict_t[period]:.3f}·C)/NH₃", lw=2
        )
        line_handles[period] = line
    ax3.set(
        title=r"$(\mathrm{NH_3}(x) - s\,\mathrm{C}(x))\,/\,\mathrm{NH_3}(x)$, $0.75<M_x^2<1.05$, $|t|<1$",
        xlabel=r"$x_{B}$",
        ylabel="ratio",
        xlim=(0.0, 1.0),
        ylim=(0.0, 1.0)
    )

    # 8) Overlay dilution factor data points with error bars, matching colors
    # Su22
    handle_su22 = line_handles["Su22"]
    color_su22 = handle_su22.get_color()
    ax3.errorbar(
        x_Su22, dil_Su22, yerr=dil_err_Su22,
        fmt='o', color=color_su22, label="Su22 dilution"
    )
    # Fa22
    handle_fa22 = line_handles["Fa22"]
    color_fa22 = handle_fa22.get_color()
    ax3.errorbar(
        x_Fa22, dil_Fa22, yerr=dil_err_Fa22,
        fmt='o', color=color_fa22, label="Fa22 dilution"
    )
    # Sp23: no data yet, so skip

    ax3.legend(loc='upper right', fontsize='small')
    plt.tight_layout()
    plt.savefig(XB_RATIO_OUTPUT)
    plt.close()


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    # 1) Load run charges (only runs ≤ MAX_RUNNUM)
    run_charges = parse_run_charges(CSV_PATH)

    # 2) Generate the 2×3 comparison plot and the 1×3 period-by-period comparison
    make_normalized_Mx2_plots(NH3_FILES, C_FILES, H2_FILES, run_charges, THREE_PANEL_OUTPUT)


if __name__ == "__main__":
    main()