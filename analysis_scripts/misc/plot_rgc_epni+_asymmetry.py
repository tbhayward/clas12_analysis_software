#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

# -----------------------------------------------------------------------------
# Hard‐coded data for all three run periods (Su22, Fa22, Sp23)
#   We only need the enpichi2FitsALUsinphi_* lists here, because the bottom‐right
#   panel extracts [row[0] for row in …] as x_B and plots the given dilution points.
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

# Dilution‐factor data you wanted to overlay in Bottom‐Right:
#   For Su22 and Fa22 we already have:
x_Su22   = np.array([row[0] for row in enpichi2FitsALUsinphi_Su22])
dil_Su22 = np.array([0.711119, 0.376873, 0.390790, 0.401177, 0.410563, 0.416077])
dil_err_Su22 = np.array([0.173083, 0.0175818, 0.00942491, 0.00801711, 0.0108668, 0.0225813])

x_Fa22   = np.array([row[0] for row in enpichi2FitsALUsinphi_Fa22])
dil_Fa22 = np.array([0.374532, 0.406101, 0.386122, 0.397153, 0.415526, 0.428239])
dil_err_Fa22 = np.array([0.117283, 0.00608866, 0.00354573, 0.00301186, 0.00400862, 0.00824607])

# (When you have Sp23 dilution data, put it in a similar form:
#   x_Sp23   = np.array([ … ])
#   dil_Sp23 = np.array([ … ])
#   dil_err_Sp23 = np.array([ … ])
# and it will automatically be plotted in green.)


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
    Returns (label, hist_no_t, hist_t, mx2_selected_t, run_base, t_vals) where:
      - hist_no_t: normalized Mx² histogram with just (runnum ≤ MAX_RUNNUM, Mx² ≤ 2.0).
      - hist_t:    normalized Mx² histogram with additional |t| < 1.0 cut.
      - mx2_selected_t: the Mx² values that survive |t|<1 (for bottom-right ratio plot).
      - run_base, t_vals: arrays of length=N_base containing run numbers and t for all base‐filtered events.
    """
    filepath, label, run_charges, quick = args

    # If QUICK_RUN=True, call the quick‐stop variant:
    if quick:
        return process_file_quick_stop(filepath, label, run_charges)

    # -----------------------------------------------------------------------------
    # Full‐array version (quick=False)
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
        # return an empty “mx2_selected_t” array for downstream logic
        return (label, hist_no_t, hist_t, np.array([]), run_base, np.array([]))

    # --- Compute “no‐t‐cut” histogram (only run/Mx² filter) ---
    run_final_no_t = run_base.copy()
    mx2_final_no_t = mx2_base.copy()

    # 5a) Select which runs contribute: all unique runs
    selected_runs_no_t = set(np.unique(run_final_no_t))

    # 6a) Sum total beam charge for those runs (Q_no_t)
    Q_no_t = 0.0
    for r in selected_runs_no_t:
        q = run_charges.get(r, 0.0)
        if q <= 0.0:
            print(f"    [WARNING] (no-t) run {r} has zero or missing charge.")
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
    run_t = run_base[mask_t]
    mx2_t = mx2_base[mask_t]

    # Count how many survive |t|<1
    num_t_events = run_t.size
    print(f"[DEBUG][{label}] after |t|<1 cut: {num_t_events} events")

    if num_t_events == 0:
        # No events survive the t-cut
        hist_t_array = np.zeros(len(bins) - 1, dtype=float)
        mx2_selected_t = np.array([])  # nothing left
    else:
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
            hist_t_array = counts_t.astype(float) / Q_t
        else:
            hist_t_array = np.zeros(len(bins) - 1, dtype=float)

        # We return the actual Mx² values that survive the t-cut, for bottom-right ratio plot
        mx2_selected_t = mx2_use_t.copy()

    # Compute survival fraction in terms of event count
    survival_frac = 100.0 * num_t_events / num_base_events if num_base_events > 0 else 0.0
    print(f"[DEBUG][{label}] event survival = {survival_frac:.1f}% ({num_t_events}/{num_base_events})")

    return (label, hist_no_t, hist_t_array, mx2_selected_t, run_base, t_vals)


def process_file_quick_stop(filepath, label, run_charges):
    """
    Quick‐stop version of process_file():
    - Only process the FIRST 5 unique runs to build both histograms.
    - Returns the same signature as process_file (though mx2_selected_t, run_base, t_vals
      will be “dummy” or partial).
    """
    bins = np.linspace(0.0, 1.5, 101)
    hist_no_t = np.zeros(len(bins) - 1, dtype=float)
    hist_t    = np.zeros(len(bins) - 1, dtype=float)

    seen_runs_no_t = set()
    seen_runs_t    = set()
    mx2_acc_no_t = []
    mx2_acc_t    = []

    total_base_events = 0
    total_t_events    = 0

    # We'll also keep a “run_base_list” and “t_vals_list” to return for completeness,
    # but they won't be the full array—just from the first few chunks.
    run_base_list = []
    t_vals_list   = []

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

        # Base filter
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

        # Loop event by event, building a small sample
        for (rnum, mx2_val, ep, eth, eph, pp, pth, pph) in zip(
                run_base_chunk, mx2_base_chunk,
                e_p_base_chunk, e_th_base_chunk, e_ph_base_chunk,
                p_p_base_chunk, p_th_base_chunk, p_ph_base_chunk
            ):

            # Collect run_base/t_vals samples
            if len(run_base_list) < 50000:  # don't store more than a handful
                run_base_list.append(rnum)
                tv = compute_t_array(
                    np.array([rnum]),
                    np.array([ep]), np.array([eth]), np.array([eph]),
                    np.array([pp]), np.array([pth]), np.array([pph])
                )[0]
                t_vals_list.append(tv)

            # “no-t” side
            if len(seen_runs_no_t) < 5:
                if rnum not in seen_runs_no_t:
                    seen_runs_no_t.add(rnum)
                if rnum in seen_runs_no_t:
                    mx2_acc_no_t.append(mx2_val)

            # “t-cut” side
            if len(seen_runs_t) < 5:
                tv = compute_t_array(
                    np.array([rnum]),
                    np.array([ep]), np.array([eth]), np.array([eph]),
                    np.array([pp]), np.array([pth]), np.array([pph])
                )[0]
                if abs(tv) < 1.0:
                    if rnum not in seen_runs_t:
                        seen_runs_t.add(rnum)
                    if rnum in seen_runs_t:
                        mx2_acc_t.append(mx2_val)
                        total_t_events += 1

            if len(seen_runs_no_t) >= 5 and len(seen_runs_t) >= 5:
                break

        if len(seen_runs_no_t) >= 5 and len(seen_runs_t) >= 5:
            break

    # Summaries
    num_base_events = total_base_events
    num_t_events    = total_t_events

    # Build “no‐t” histogram
    seen_runs_list_no_t = list(seen_runs_no_t)
    Q_no_t = sum(run_charges.get(r, 0.0) for r in seen_runs_list_no_t)
    if Q_no_t > 0.0 and len(mx2_acc_no_t) > 0:
        counts_no_t, _ = np.histogram(np.array(mx2_acc_no_t), bins=bins)
        hist_no_t      = counts_no_t.astype(float) / Q_no_t

    # Build “t‐cut” histogram
    seen_runs_list_t = list(seen_runs_t)
    Q_t = sum(run_charges.get(r, 0.0) for r in seen_runs_list_t)
    if Q_t > 0.0 and len(mx2_acc_t) > 0:
        counts_t, _ = np.histogram(np.array(mx2_acc_t), bins=bins)
        hist_t      = counts_t.astype(float) / Q_t

    # Debug print
    survival_frac = 100.0 * num_t_events / num_base_events if num_base_events > 0 else 0.0
    print(f"[DEBUG][{label}, quick] base_e = {num_base_events}  |  t_e = {num_t_events}  |  "
          f"survival = {survival_frac:.1f}% ({num_t_events}/{num_base_events})")

    # Return a partial “run_base” and “t_vals” as 1D arrays
    return (label, hist_no_t, hist_t, np.array(mx2_acc_t), np.array(run_base_list), np.array(t_vals_list))


# -----------------------------------------------------------------------------
# PLOTTING FUNCTION
# -----------------------------------------------------------------------------

def make_normalized_Mx2_plots(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 2×2 figure:

      Top‐Left   (0,0): NH₃ + H₂   (no |t|-cut)
      Top‐Right  (0,1): C   + H₂   (no |t|-cut)
      Bot‐Left   (1,0): (NH₃ – C) + H₂   (with |t|<1)
      Bot‐Right  (1,1): [ (NH₃ – C)/(NH₃) ] vs x_B
                       for 0.75 < Mx² < 1.05  & |t|<1

    Uses ProcessPoolExecutor to parallelize histogramming of each file.
    """
    # 1) Prepare tasks
    tasks = []
    for fp, lbl in nh3_files + c_files + h2_files:
        tasks.append((fp, lbl, run_charges, QUICK_RUN))

    # 2) Submit tasks in parallel
    results = {}  # label → (hist_no_t, hist_t, mx2_selected_t, run_base, t_vals)
    with ProcessPoolExecutor() as executor:
        future_to_label = {executor.submit(process_file, args): args[1] for args in tasks}
        for future in as_completed(future_to_label):
            lbl = future_to_label[future]
            try:
                label, hist_no_t, hist_t, mx2_selected_t, run_base, t_vals = future.result()
                results[label] = (hist_no_t, hist_t, mx2_selected_t, run_base, t_vals)
            except Exception as e:
                print(f"[ERROR] processing {lbl}: {e}")
                # fallback zero‐arrays
                results[lbl] = (np.zeros(100,), np.zeros(100,), np.array([]), np.array([]), np.array([]))

    # 3) Prepare bins & bin centers
    bins = np.linspace(0.0, 1.5, 101)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    # 4) Create the 2×2 figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharey=False)

    ## ────────────────────────────────────────────────────────────────────────────
    # Top‐Left (0,0): NH₃ + H₂ (no |t|-cut)
    ## ────────────────────────────────────────────────────────────────────────────
    top_left_vals = []
    for _, lbl in nh3_files:
        hist_no_t, _, _, _, _ = results.get(lbl, (np.zeros(100,),)*5)
        axes[0, 0].step(bins[:-1], hist_no_t, where='post', label=lbl)
        top_left_vals.append(hist_no_t.max())

    h2_label = h2_files[0][1]
    hist_h2_no_t, _, _, _, _ = results.get(h2_label, (np.zeros(100,),)*5)
    axes[0, 0].step(
        bins[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    top_left_vals.append(hist_h2_no_t.max())

    axes[0, 0].set(
        title="NH₃ + H₂: Normalized Mx² (no |t|‐cut)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    axes[0, 0].legend(loc='upper right', fontsize='small')

    y_top_left = 1.1 * max(top_left_vals) if top_left_vals else 0.01
    axes[0, 0].set_ylim(0.0, y_top_left)

    ## ────────────────────────────────────────────────────────────────────────────
    # Top‐Right (0,1): C + H₂ (no |t|-cut)
    ## ────────────────────────────────────────────────────────────────────────────
    top_right_vals = []
    for _, lbl in c_files:
        hist_no_t, _, _, _, _ = results.get(lbl, (np.zeros(100,),)*5)
        axes[0, 1].step(bins[:-1], hist_no_t, where='post', label=lbl)
        top_right_vals.append(hist_no_t.max())

    axes[0, 1].step(
        bins[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    top_right_vals.append(hist_h2_no_t.max())

    axes[0, 1].set(
        title="C + H₂: Normalized Mx² (no |t|‐cut)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    axes[0, 1].legend(loc='upper right', fontsize='small')

    y_top_right = 1.1 * max(top_right_vals) if top_right_vals else 0.01
    axes[0, 1].set_ylim(0.0, y_top_right)

    ## ────────────────────────────────────────────────────────────────────────────
    # Bottom‐Left (1,0): (NH₃ – C) + H₂ (|t|<1)
    ## ────────────────────────────────────────────────────────────────────────────
    bottom_left_vals = []
    # Build period‐level histograms from the “with‐t” results:
    period_names = ["Su22", "Fa22", "Sp23"]

    # Collect NH₃ with‐t and C with‐t histograms:
    nh3_hist_t = {lbl: results.get(lbl, (None, np.zeros(100,), None, None, None))[1]
                  for _, lbl in nh3_files}
    c_hist_t   = {lbl: results.get(lbl, (None, np.zeros(100,), None, None, None))[1]
                  for _, lbl in c_files}

    # Compute scaled difference NH₃ – C (|t|<1) for each period:
    diff_t = {}
    for period_label in period_names:
        nh3_lbl = f"{period_label}-NH3"
        c_lbl   = f"{period_label}-C"
        nh3_vals = nh3_hist_t.get(nh3_lbl, np.zeros(100,))
        c_vals   = c_hist_t.get(c_lbl,   np.zeros(100,))

        mask_0_5 = (bin_centers >= 0.0) & (bin_centers < 0.5)
        nh3_sum = nh3_vals[mask_0_5].sum()
        c_sum   = c_vals[mask_0_5].sum()
        scale = (nh3_sum / c_sum) if c_sum > 0.0 else 0.0

        print(f"[DEBUG] (|t|<1) {period_label}: NH₃ sum(0–0.5)={nh3_sum:.3f}, "
              f"C sum(0–0.5)={c_sum:.3f}, scale={scale:.3f}")

        diff_t[period_label] = nh3_vals - (c_vals * scale)

    for idx, period_label in enumerate(period_names):
        diff = diff_t.get(period_label, np.zeros(100,))
        style = ['-', '--', '-.', ':'][idx % 4]
        axes[1, 0].step(
            bins[:-1], diff, where='post',
            linestyle=style, label=f"{period_label} NH₃–C"
        )
        bottom_left_vals.append(diff.max())

    # Also plot H₂ with‐t
    _, hist_h2_t, _, _, _ = results.get(h2_label, (None, np.zeros(100,), None, None, None))
    axes[1, 0].step(
        bins[:-1], hist_h2_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    bottom_left_vals.append(hist_h2_t.max())

    axes[1, 0].set(
        title="(NH₃ – C) + H₂ (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    axes[1, 0].legend(loc='upper right', fontsize='small')

    y_bottom_left = 1.1 * max(bottom_left_vals) if bottom_left_vals else 0.01
    axes[1, 0].set_ylim(0.0, y_bottom_left)

    ## ────────────────────────────────────────────────────────────────────────────
    # Bottom‐Right (1,1): [ (NH₃–C)/(NH₃) ] vs x_B  for 0.75<Mx²<1.05 & |t|<1
    ## ────────────────────────────────────────────────────────────────────────────

    # We need to reconstruct, for each period, a small “normalized yield” histogram
    # on a fine grid of x_B.  But in practice we’ll do it “point‐by‐point” by scanning
    # all events that survived the t‐cut, selecting only those with 0.75 < Mx² < 1.05,
    # then computing (NH₃ – C)/NH₃ in small bins of x_B.

    # Step A: gather all “MX² & t‐cut survivors” + their x_B from each period’s NH3 and C
    # We do that by re‐opening each ROOT file and selecting:
    #    (run ≤ MAX_RUNNUM,  0.75<Mx²<1.05,  |t|<1), then computing t & x_B per event.

    def compute_xB_array(run_arr, e_p_arr, e_th_arr, e_ph_arr):
        """
        Given arrays (run, e_p, e_th, e_ph) for exclusive ep→e'πn,
        compute x_B ≡ Q²/(2M_pν), where ν = E_beam – E_e'.
        Q² = 4 E_beam E_e' sin²(θ_e/2).
        (This formula holds for electron scattering, with masses in GeV.)
        """
        Eb = get_beam_energy(run_arr)
        E_e = np.sqrt(e_p_arr**2 + m_e**2)
        ν   = Eb - E_e
        sin2 = np.sin(e_th_arr / 2.0)**2
        Q2  = 4 * Eb * E_e * sin2
        # Avoid division by zero or negative ν:
        mask_valid = (ν > 0)
        xB = np.zeros_like(run_arr, dtype=float)
        xB[mask_valid] = Q2[mask_valid] / (2.0 * m_p * ν[mask_valid])
        return xB

    # For each period ("Su22", "Fa22", "Sp23"), we will:
    #  1) Re-open NH₃ TTree and C TTree
    #  2) Apply the “run≤MAX_RUNNUM & 0.75<Mx²<1.05” mask
    #  3) Compute t and x_B for each surviving event
    #  4) Keep only |t|<1
    #  5) Now we have two arrays of x_B: xB_NH3 and xB_C
    #  6) Bin those x_B into, say, 25 equally spaced bins from x=0→0.6 (you can adjust)
    #  7) Form histograms H_N(x) and H_C(x), then compute “(H_N – H_C)/H_N” bin by bin.

    # We’ll define a helper function to do exactly that:
    def get_ratio_NminusC_over_N(period_label, mh2_allowed=False):
        """
        For a given period_label like "Su22", return (bin_centers, ratio, ratio_err)
        where ratio = (H_NH3(x) – H_C(x))/H_NH3(x), for events that satisfy:
          * run ≤ MAX_RUNNUM
          * 0.75 < Mx² < 1.05
          * |t| < 1
        Binnings: 25 bins between x_B=0.0 and x_B=0.6

        Returns:
          x_centers (length=25), ratio_array (length=25), ratio_err_array (length=25)
        """
        # 1) pick filepaths:
        fp_nh3 = next(fp for fp, lbl in nh3_files if lbl.startswith(period_label))
        fp_c   = next(fp for fp, lbl in c_files   if lbl.startswith(period_label))

        # 2) read in arrays chunk‐by‐chunk for NH₃:
        x_nh3_list = []
        x_c_list   = []

        # We’ll do both in two separate uproot.iterate loops (they’re large).
        # (a) NH₃:
        for chunk in uproot.iterate(
            fp_nh3 + ":PhysicsEvents",
            ["runnum", "Mx2", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi"],
            step_size=200_000
        ):
            run_arr  = np.array(chunk["runnum"], dtype=int)
            mx2_arr  = np.array(chunk["Mx2"], dtype=float)
            e_p_arr  = np.array(chunk["e_p"], dtype=float)
            e_th_arr = np.array(chunk["e_theta"], dtype=float)
            e_ph_arr = np.array(chunk["e_phi"], dtype=float)
            p_p_arr  = np.array(chunk["p_p"], dtype=float)
            p_th_arr = np.array(chunk["p_theta"], dtype=float)
            p_ph_arr = np.array(chunk["p_phi"], dtype=float)

            # mask: run ≤ MAX_RUNNUM, 0.75 < Mx² < 1.05
            mask_period = (
                (run_arr <= MAX_RUNNUM) &
                (mx2_arr >  0.75) &
                (mx2_arr <  1.05)
            )
            run_base_p  = run_arr[mask_period]
            e_p_p       = e_p_arr[mask_period]
            e_th_p      = e_th_arr[mask_period]
            e_ph_p      = e_ph_arr[mask_period]
            p_p_p       = p_p_arr[mask_period]
            p_th_p      = p_th_arr[mask_period]
            p_ph_p      = p_ph_arr[mask_period]

            if run_base_p.size == 0:
                continue

            # compute t for these events:
            tvals_p = compute_t_array(
                run_base_p,
                e_p_p, e_th_p, e_ph_p,
                p_p_p, p_th_p, p_ph_p
            )
            mask_tcut = np.abs(tvals_p) < 1.0

            if not np.any(mask_tcut):
                continue

            # compute x_B for events that survive t-cut
            run_tcut = run_base_p[mask_tcut]
            e_p_tcut  = e_p_p[mask_tcut]
            e_th_tcut = e_th_p[mask_tcut]
            e_ph_tcut = e_ph_p[mask_tcut]

            xB_vals = compute_xB_array(run_tcut, e_p_tcut, e_th_tcut, e_ph_tcut)
            x_nh3_list.append(xB_vals)

        # (b) C:
        for chunk in uproot.iterate(
            fp_c + ":PhysicsEvents",
            ["runnum", "Mx2", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi"],
            step_size=200_000
        ):
            run_arr  = np.array(chunk["runnum"], dtype=int)
            mx2_arr  = np.array(chunk["Mx2"], dtype=float)
            e_p_arr  = np.array(chunk["e_p"], dtype=float)
            e_th_arr = np.array(chunk["e_theta"], dtype=float)
            e_ph_arr = np.array(chunk["e_phi"], dtype=float)
            p_p_arr  = np.array(chunk["p_p"], dtype=float)
            p_th_arr = np.array(chunk["p_theta"], dtype=float)
            p_ph_arr = np.array(chunk["p_phi"], dtype=float)

            mask_period = (
                (run_arr <= MAX_RUNNUM) &
                (mx2_arr >  0.75) &
                (mx2_arr <  1.05)
            )
            run_base_p  = run_arr[mask_period]
            e_p_p       = e_p_arr[mask_period]
            e_th_p      = e_th_arr[mask_period]
            e_ph_p      = e_ph_arr[mask_period]
            p_p_p       = p_p_arr[mask_period]
            p_th_p      = p_th_arr[mask_period]
            p_ph_p      = p_ph_arr[mask_period]

            if run_base_p.size == 0:
                continue

            tvals_p = compute_t_array(
                run_base_p,
                e_p_p, e_th_p, e_ph_p,
                p_p_p, p_th_p, p_ph_p
            )
            mask_tcut = np.abs(tvals_p) < 1.0

            if not np.any(mask_tcut):
                continue

            e_p_tcut  = e_p_p[mask_tcut]
            e_th_tcut = e_th_p[mask_tcut]
            e_ph_tcut = e_ph_p[mask_tcut]
            run_tcut  = run_base_p[mask_tcut]

            xB_vals = compute_xB_array(run_tcut, e_p_tcut, e_th_tcut, e_ph_tcut)
            x_c_list.append(xB_vals)

        # Concatenate all chunks
        if len(x_nh3_list) == 0:
            x_nh3 = np.array([])
        else:
            x_nh3 = np.concatenate(x_nh3_list)

        if len(x_c_list) == 0:
            x_c = np.array([])
        else:
            x_c = np.concatenate(x_c_list)

        # If there are no events in either NH₃ or C, return empties
        if x_nh3.size == 0 or x_c.size == 0:
            # We'll return 25 bins of zeros
            x_centers = 0.5 * (np.linspace(0.0, 0.6, 26)[:-1] + np.linspace(0.0, 0.6, 26)[1:])
            return (x_centers, np.zeros_like(x_centers), np.zeros_like(x_centers))

        # 6) Bin these x_B arrays into 25 bins from 0→0.6:
        x_bins = np.linspace(0.0, 0.6, 26)  # 25 bins
        H_n, _ = np.histogram(x_nh3, bins=x_bins)
        H_c, _ = np.histogram(x_c,   bins=x_bins)

        # Convert to “events/nC” shape—BUT since we only care about ratio (H_n - H_c)/H_n,
        # the absolute normalization cancels. We just need to compute bin‐by‐bin ratio.
        # We also want an approximate “stat‐error” on that ratio: if H_n is large, the 
        # binomial error on C is sqrt(H_c), etc.  For simplicity, we’ll estimate:
        #   ratio_i = (H_n[i] - H_c[i]) / H_n[i]
        #   error_i = sqrt( (sqrt(H_c[i]))² + ( (H_c[i]/H_n[i]) * sqrt(H_n[i]) )² ) / H_n[i]
        # (propagate errors of H_n and H_c → ratio). If H_n=0, skip that bin.

        x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
        ratio    = np.zeros_like(x_centers)
        ratio_e  = np.zeros_like(x_centers)

        for i in range(len(x_centers)):
            n_i = H_n[i]
            c_i = H_c[i]
            if n_i <= 0:
                ratio[i] = 0.0
                ratio_e[i] = 0.0
            else:
                val = (n_i - c_i) / n_i
                # error: 
                #   σ(H_n)=√n_i,  σ(H_c)=√c_i
                #   using formula: f = (n - c)/n = 1 - (c/n)
                #   var(f) ≈ var(c)/n² + (c² var(n))/n⁴
                #   = (c)/n² + c²/(n⁴) * n
                #   = (c)/n² + c²/(n³)
                # do sqrt:
                term1 = c_i / (n_i**2)
                term2 = (c_i**2) / (n_i**3)
                ratio[i]   = val
                ratio_e[i] = np.sqrt(term1 + term2)
        return (x_centers, ratio, ratio_e)

    # Now call that helper for each of the three periods:
    x_su22,  r_su22,  e_su22  = get_ratio_NminusC_over_N("Su22")
    x_fa22,  r_fa22,  e_fa22  = get_ratio_NminusC_over_N("Fa22")
    x_sp23,  r_sp23,  e_sp23  = get_ratio_NminusC_over_N("Sp23")

    # Finally, plot all four panels:

    # -------------------------------------------------------------------------
    # Panel (0,0): NH₃ + H₂ (no |t|‐cut)
    # -------------------------------------------------------------------------
    # Already done above.

    # -------------------------------------------------------------------------
    # Panel (0,1): C + H₂ (no |t|‐cut)
    # -------------------------------------------------------------------------
    # Already done above.

    # -------------------------------------------------------------------------
    # Panel (1,0): (NH₃ – C) + H₂ (|t|<1)
    # -------------------------------------------------------------------------
    # Already done above.

    # -------------------------------------------------------------------------
    # Panel (1,1): [ (NH₃ – C)/(NH₃) ] vs x_B
    # -------------------------------------------------------------------------
    # Plot each period’s curve:
    axes[1, 1].plot(x_su22, r_su22,   '-', color='C0', label="Su22")
    axes[1, 1].plot(x_fa22, r_fa22,  '--', color='C1', label="Fa22")
    axes[1, 1].plot(x_sp23, r_sp23, '-.', color='C2', label="Sp23")

    # Overlay Su22 dilution points (blue) and Fa22 (orange):
    axes[1, 1].errorbar(
        x_Su22, dil_Su22, yerr=dil_err_Su22,
        fmt='o', color='C0', capsize=3, label="Su22 dil. data"
    )
    axes[1, 1].errorbar(
        x_Fa22, dil_Fa22, yerr=dil_err_Fa22,
        fmt='s', color='C1', capsize=3, label="Fa22 dil. data"
    )
    # (If you have x_Sp23, dil_Sp23, dil_err_Sp23, add a similar:
    #   axes[1,1].errorbar(x_Sp23, dil_Sp23, yerr=dil_err_Sp23,
    #                      fmt='^', color='C2', capsize=3, label="Sp23 dil. data")
    # )

    axes[1, 1].set(
        title=r"$(\mathrm{NH_3}-\mathrm{C})/\mathrm{NH_3}$ vs $x_B$",
        xlabel=r"$x_B$",
        ylabel=r"$(N - C)/N$",
        xlim=(0.0, 0.6)
    )
    axes[1, 1].legend(loc='upper left', fontsize='small')

    # Adjust y‐limits so that the data points fit nicely
    all_vals = np.concatenate([r_su22, r_fa22])
    if all_vals.size > 0:
        y_hi = 1.1 * np.nanmax(all_vals)
        axes[1, 1].set_ylim(0.0, y_hi)

    # -------------------------------------------------------------------------
    # Finally, save the 2×2 figure
    # -------------------------------------------------------------------------
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    run_charges = parse_run_charges(CSV_PATH)
    make_normalized_Mx2_plots(NH3_FILES, C_FILES, H2_FILES, run_charges, THREE_PANEL_OUTPUT)


if __name__ == "__main__":
    main()