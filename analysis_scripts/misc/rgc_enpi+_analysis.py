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
OUTPUT_DIR         = "output/enpi+"
FOUR_PANEL_OUTPUT  = os.path.join(OUTPUT_DIR, "four_panel_Mx2_x_comparison.pdf")

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
    Read clas12_run_info.csv and return a dict { runnum: charge }, skipping any run > MAX_RUNNUM.
    Lines beginning with '#' are ignored.
    """
    run_charges = {}
    if not os.path.exists(csv_path):
        print(f"[ERROR] CSV not found: {csv_path}")
        return run_charges

    with open(csv_path, "r") as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith('#'):
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


# -----------------------------------------------------------------------------
# PROCESSING FUNCTIONS
# -----------------------------------------------------------------------------

def process_file(args):
    """
    Worker for parallel histogramming of Mx²:
    Returns (label, hist_no_t, hist_t) where:
      - hist_no_t: normalized Mx² histogram with just (run ≤ MAX_RUNNUM, Mx² ≤ 2.0).
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

    # --- Compute “no-t-cut” histogram (run/Mx² filter only) ---
    run_final_no_t = run_base.copy()
    mx2_final_no_t = mx2_base.copy()

    # 5a) Which runs contribute (quick=False ⇒ all unique)
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

    # --- Compute “t-cut” histogram (also require |t|<1) ---
    t_vals = compute_t_array(
        run_base,
        e_p_base, e_th_base, e_ph_base,
        p_p_base, p_th_base, p_ph_base
    )

    mask_t = np.abs(t_vals) < 1.0
    run_t  = run_base[mask_t]
    mx2_t  = mx2_base[mask_t]

    num_t_events = run_t.size
    print(f"[DEBUG][{label}] after |t|<1 cut: {num_t_events} events")

    if num_t_events > 0:
        # 5b) Which runs survive under |t|<1
        selected_runs_t = set(np.unique(run_t))

        # 6b) Sum total beam charge for those runs (Q_t)
        Q_t = 0.0
        for r in selected_runs_t:
            q = run_charges.get(r, 0.0)
            if q <= 0.0:
                print(f"    [WARNING] (t) run {r} has zero or missing charge.")
            Q_t += q

        # 7b) Mask t-filtered events so that only those runs remain
        mask_runs_t = np.isin(run_t, list(selected_runs_t))
        mx2_use_t   = mx2_t[mask_runs_t]

        if Q_t > 0.0 and mx2_use_t.size > 0:
            counts_t, _ = np.histogram(mx2_use_t, bins=bins)
            hist_t      = counts_t.astype(float) / Q_t

    # Compute survival fraction in terms of raw event count
    survival_frac = 100.0 * num_t_events / num_base_events if num_base_events > 0 else 0.0
    print(f"[DEBUG][{label}] event survival = {survival_frac:.1f}% ({num_t_events}/{num_base_events})")

    return (label, hist_no_t, hist_t)


def process_file_quick_stop(filepath, label, run_charges):
    """
    “Quick-stop” version of process_file():
    - Only process the FIRST 5 unique runs to build both histograms.
    - Uses uproot.iterate() to read in chunks, stops once 5 distinct runs appear in each category.
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
    total_t_events    = 0

    # Iterate in chunks of 100 k events, stop once each set has 5 distinct runs
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

        # Now loop event-by-event until 5 distinct runs appear in both categories
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

            # (2) Compute t for this event if we still need runs for the t-cut
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

            # If both categories already have 5 runs, break
            if (len(seen_runs_no_t) >= 5) and (len(seen_runs_t) >= 5):
                break

        if (len(seen_runs_no_t) >= 5) and (len(seen_runs_t) >= 5):
            break

    # We know total_base_events (though we only recorded Mx² from—at most—the first 5 runs),
    # and mx2_acc_no_t, mx2_acc_t hold Mx² from “those first 5 runs” in each category.

    num_base_events = total_base_events
    num_no_t_events = len(mx2_acc_no_t)
    num_t_events    = len(mx2_acc_t)

    # Build “no-t” histogram
    seen_runs_list_no_t = list(seen_runs_no_t)
    Q_no_t = sum(run_charges.get(r, 0.0) for r in seen_runs_list_no_t)
    if Q_no_t > 0.0 and num_no_t_events > 0:
        counts_no_t, _ = np.histogram(np.array(mx2_acc_no_t), bins=bins)
        hist_no_t      = counts_no_t.astype(float) / Q_no_t
    else:
        hist_no_t = np.zeros(len(bins) - 1, dtype=float)

    # Build “t-cut” histogram
    seen_runs_list_t = list(seen_runs_t)
    Q_t = sum(run_charges.get(r, 0.0) for r in seen_runs_list_t)
    if Q_t > 0.0 and num_t_events > 0:
        counts_t, _ = np.histogram(np.array(mx2_acc_t), bins=bins)
        hist_t      = counts_t.astype(float) / Q_t
    else:
        hist_t = np.zeros(len(bins) - 1, dtype=float)

    survival_frac = 100.0 * num_t_events / num_base_events if num_base_events > 0 else 0.0
    print(f"[DEBUG][{label}, quick] base_e = {num_base_events}  |  t_e = {num_t_events}  |  "
          f"survival = {survival_frac:.1f}% ({num_t_events}/{num_base_events})")

    return (label, hist_no_t, hist_t)


# -----------------------------------------------------------------------------
# PLOTTING FUNCTION: 2×2 canvas
# -----------------------------------------------------------------------------

def make_four_panel_plots(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 2×2 figure (no shared y):
      ┌─────────────────────────────────────────┐ ┌─────────────────────────────────────────┐
      │ (0,0) NH₃ + H₂: Mx² (no t-cut)         │ │ (0,1)  C + H₂: Mx² (no t-cut)          │
      │               (step(hist_no_t))        │ │               (step(hist_no_t))        │
      └─────────────────────────────────────────┘ └─────────────────────────────────────────┘
      ┌─────────────────────────────────────────┐ ┌─────────────────────────────────────────┐
      │ (1,0) [NH₃–C] + H₂: Mx² (no t-cut)     │ │ (1,1) [NH₃(x)−C(x)] / NH₃(x)            │
      │               (step(diff_no_t))        │ │   0.75<Mx²<1.05, |t|<1, bin in x        │
      └─────────────────────────────────────────┘ └─────────────────────────────────────────┘

    Uses ProcessPoolExecutor to parallelize the first three Mx² histograms,
    then does a separate x–distribution for panel (1,1).
    """
    # 1) Build tasks for parallel Mx² histogramming
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
                # On error, store zero arrays
                results[lbl] = (np.zeros(100, dtype=float), np.zeros(100, dtype=float))

    # 3) Prepare bins & bin centers for Mx²
    bins_mx2     = np.linspace(0.0, 1.5, 101)
    centers_mx2  = 0.5 * (bins_mx2[:-1] + bins_mx2[1:])

    # 4) Start the 2×2 figure (no shared y)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharey=False)
    axes = axes.flatten()

    # Keep track of maxima to set dynamic y‐limits:
    top_row_vals    = []
    bottom_left_vals = []

    ## ────────────────────────────────────────────────────────────────────────────
    # Panel (0,0): NH₃ + H₂ (no t-cut)
    ## ────────────────────────────────────────────────────────────────────────────
    ax00 = axes[0]
    for _, lbl in nh3_files:
        hist_no_t, _ = results.get(lbl, (np.zeros(100), np.zeros(100)))
        ax00.step(bins_mx2[:-1], hist_no_t, where='post', label=lbl)
        top_row_vals.append(hist_no_t.max())

    h2_label = h2_files[0][1]
    hist_h2_no_t, _ = results.get(h2_label, (np.zeros(100), np.zeros(100)))
    ax00.step(
        bins_mx2[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    top_row_vals.append(hist_h2_no_t.max())

    ax00.set(
        title="NH₃ + H₂: Normalized $M_{x}^2$ (no $t$-cut)",
        xlabel=r"$M_{x}^2$ (GeV$^2$)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    ax00.legend(loc='upper right', fontsize='small')

    ## ────────────────────────────────────────────────────────────────────────────
    # Panel (0,1): C + H₂ (no t-cut)
    ## ────────────────────────────────────────────────────────────────────────────
    ax01 = axes[1]
    for _, lbl in c_files:
        hist_no_t, _ = results.get(lbl, (np.zeros(100), np.zeros(100)))
        ax01.step(bins_mx2[:-1], hist_no_t, where='post', label=lbl)
        top_row_vals.append(hist_no_t.max())

    ax01.step(
        bins_mx2[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    top_row_vals.append(hist_h2_no_t.max())

    ax01.set(
        title="C + H₂: Normalized $M_{x}^2$ (no $t$-cut)",
        xlabel=r"$M_{x}^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    ax01.legend(loc='upper right', fontsize='small')

    # Determine a common y‐limit for the entire top row:
    y_top = 1.1 * max(top_row_vals) if top_row_vals else 0.1
    ax00.set(ylim=(0.0, y_top))
    ax01.set(ylim=(0.0, y_top))

    ## ────────────────────────────────────────────────────────────────────────────
    # Panel (1,0): [NH₃ – C] + H₂ (no t-cut)
    ## ────────────────────────────────────────────────────────────────────────────
    ax10 = axes[2]
    period_names = ["Su22", "Fa22", "Sp23"]

    # Gather NH₃ and C no-t histograms into dictionaries
    nh3_hist_no_t = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[0] for _, lbl in nh3_files}
    c_hist_no_t   = {lbl: results.get(lbl, (np.zeros(100), np.zeros(100)))[0] for _, lbl in c_files}

    diff_no_t = {}
    for period_label in period_names:
        nh3_lbl = f"{period_label}-NH3"
        c_lbl   = f"{period_label}-C"
        nh3_vals = nh3_hist_no_t.get(nh3_lbl, np.zeros(100))
        c_vals   = c_hist_no_t.get(c_lbl,   np.zeros(100))

        # Scale C so that NH₃ and C match in 0 ≤ Mx² < 0.5
        mask_0_5 = (centers_mx2 >= 0.0) & (centers_mx2 < 0.5)
        nh3_sum = nh3_vals[mask_0_5].sum()
        c_sum   = c_vals[mask_0_5].sum()
        scale = (nh3_sum / c_sum) if c_sum > 0.0 else 0.0
        print(f"[DEBUG][no-t] {period_label}: NH₃ sum(0–0.5)={nh3_sum:.3f}, C sum(0–0.5)={c_sum:.3f}, scale={scale:.3f}")
        diff_no_t[period_label] = nh3_vals - (c_vals * scale)

    for idx, period_label in enumerate(period_names):
        diff = diff_no_t.get(period_label, np.zeros(100))
        style = ['-', '--', '-.', ':'][idx % 4]
        ax10.step(
            bins_mx2[:-1], diff, where='post',
            linestyle=style, label=f"{period_label} NH₃–C"
        )
        bottom_left_vals.append(diff.max())

    # Overplot H₂ (no-t) on the same panel
    ax10.step(
        bins_mx2[:-1], hist_h2_no_t, where='post',
        color='k', linestyle='-', label=h2_label
    )
    bottom_left_vals.append(hist_h2_no_t.max())

    ax10.set(
        title="(NH₃ – C) + H₂ (no $t$-cut)",
        xlabel=r"$M_{x}^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    ax10.legend(loc='upper right', fontsize='small')

    ## ────────────────────────────────────────────────────────────────────────────
    # Panel (1,1): Ratio versus xB under 0.75<Mx²<1.05, |t|<1
    ## ────────────────────────────────────────────────────────────────────────────
    ax11 = axes[3]

    # We will build NH₃(x) and C(x) histograms in the same run, under the simultaneous cuts:
    #   0.75 < Mx² < 1.05,   |t| < 1.0
    # Then plot [NH₃(x) - C(x)]/NH₃(x)  vs  x.
    #
    # We'll use 50 uniform bins in x ∈ [0,1].  (You may adjust as needed.)
    x_bins = np.linspace(0.0, 1.0, 51)
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

    # Initialize “accumulated” arrays for NH₃ and C
    #   – We want counts per x-bin, then normalize by total charge so that NH₃(x), C(x) 
    #     each become “events / nC”.
    nh3_counts_x = np.zeros(len(x_bins) - 1, dtype=float)
    c_counts_x   = np.zeros(len(x_bins) - 1, dtype=float)

    nh3_charge_total = 0.0
    c_charge_total   = 0.0

    # ---- Build NH₃(x) under 0.75<Mx²<1.05, |t|<1 ----
    for filepath, label in nh3_files:
        tree = uproot.open(filepath)["PhysicsEvents"]

        # Extract all arrays in one go
        run_arr  = tree["runnum"].array(library="np").astype(int)
        mx2_arr  = tree["Mx2"].array(library="np")
        x_arr    = tree["x"].array(library="np")
        e_p_arr  = tree["e_p"].array(library="np")
        e_th_arr = tree["e_theta"].array(library="np")
        e_ph_arr = tree["e_phi"].array(library="np")
        p_p_arr  = tree["p_p"].array(library="np")
        p_th_arr = tree["p_theta"].array(library="np")
        p_ph_arr = tree["p_phi"].array(library="np")

        # 1) Baseline filter: run ≤ MAX_RUNNUM AND 0.75 < Mx² < 1.05
        mask_base_x = (run_arr <= MAX_RUNNUM) & (mx2_arr > 0.75) & (mx2_arr < 1.05)
        run_base_x  = run_arr[mask_base_x]
        mx2_base_x  = mx2_arr[mask_base_x]
        x_base_x    = x_arr[mask_base_x]
        e_p_base_x  = e_p_arr[mask_base_x]
        e_th_base_x = e_th_arr[mask_base_x]
        e_ph_base_x = e_ph_arr[mask_base_x]
        p_p_base_x  = p_p_arr[mask_base_x]
        p_th_base_x = p_th_arr[mask_base_x]
        p_ph_base_x = p_ph_arr[mask_base_x]

        if run_base_x.size == 0:
            continue

        # 2) Compute t-array for these base-filtered events
        t_vals_x = compute_t_array(
            run_base_x,
            e_p_base_x, e_th_base_x, e_ph_base_x,
            p_p_base_x, p_th_base_x, p_ph_base_x
        )
        mask_t_x = np.abs(t_vals_x) < 1.0

        run_t_x = run_base_x[mask_t_x]
        x_t_x   = x_base_x[mask_t_x]

        if run_t_x.size == 0:
            continue

        # 3) Among those |t|<1 events, find the unique runs that contribute
        unique_runs = np.unique(run_t_x)
        for r in unique_runs:
            q = run_charges.get(r, 0.0)
            if q <= 0.0:
                print(f"    [WARNING][NH₃,x] run {r} has zero or missing charge.")
            nh3_charge_total += q

        # 4) Now fill “x” histogram *only* for events whose run ∈ unique_runs
        mask_run_select = np.isin(run_t_x, unique_runs)
        x_select = x_t_x[mask_run_select]
        if x_select.size > 0:
            counts_x, _ = np.histogram(x_select, bins=x_bins)
            nh3_counts_x += counts_x

    # ---- Build C(x) under 0.75<Mx²<1.05, |t|<1 ----
    for filepath, label in c_files:
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

        mask_base_x = (run_arr <= MAX_RUNNUM) & (mx2_arr > 0.75) & (mx2_arr < 1.05)
        run_base_x  = run_arr[mask_base_x]
        mx2_base_x  = mx2_arr[mask_base_x]
        x_base_x    = x_arr[mask_base_x]
        e_p_base_x  = e_p_arr[mask_base_x]
        e_th_base_x = e_th_arr[mask_base_x]
        e_ph_base_x = e_ph_arr[mask_base_x]
        p_p_base_x  = p_p_arr[mask_base_x]
        p_th_base_x = p_th_arr[mask_base_x]
        p_ph_base_x = p_ph_arr[mask_base_x]

        if run_base_x.size == 0:
            continue

        t_vals_x = compute_t_array(
            run_base_x,
            e_p_base_x, e_th_base_x, e_ph_base_x,
            p_p_base_x, p_th_base_x, p_ph_base_x
        )
        mask_t_x = np.abs(t_vals_x) < 1.0

        run_t_x = run_base_x[mask_t_x]
        x_t_x   = x_base_x[mask_t_x]

        if run_t_x.size == 0:
            continue

        unique_runs = np.unique(run_t_x)
        for r in unique_runs:
            q = run_charges.get(r, 0.0)
            if q <= 0.0:
                print(f"    [WARNING][C,x] run {r} has zero or missing charge.")
            c_charge_total += q

        mask_run_select = np.isin(run_t_x, unique_runs)
        x_select = x_t_x[mask_run_select]
        if x_select.size > 0:
            counts_x, _ = np.histogram(x_select, bins=x_bins)
            c_counts_x += counts_x

    # 5) Now build “normalized” NH₃(x) and C(x):
    nh3_norm_x = np.zeros_like(nh3_counts_x, dtype=float)
    c_norm_x   = np.zeros_like(c_counts_x, dtype=float)

    if nh3_charge_total > 0.0:
        nh3_norm_x = nh3_counts_x.astype(float) / nh3_charge_total
    if c_charge_total > 0.0:
        c_norm_x = c_counts_x.astype(float) / c_charge_total

    # 6) Finally compute ratio:  [NH₃(x) – C(x)] / NH₃(x)  bin-by-bin
    ratio_x = np.zeros_like(nh3_norm_x, dtype=float)
    for i in range(len(ratio_x)):
        if nh3_norm_x[i] > 0.0:
            ratio_x[i] = (nh3_norm_x[i] - c_norm_x[i]) / nh3_norm_x[i]
        else:
            ratio_x[i] = 0.0

    # Plot ratio on ax11
    ax11.step(x_centers, ratio_x, where='mid', color='m', label=r"$(\mathrm{NH_3}-\mathrm{C})/\mathrm{NH_3}$")
    ax11.set(
        title=r"$(\mathrm{NH_3}(x)-\mathrm{C}(x))\,/\,\mathrm{NH_3}(x)$, $0.75<M_x^2<1.05$, $|t|<1$",
        xlabel=r"$x_{B}$",
        ylabel="ratio",
        xlim=(0.0, 1.0),
        ylim=(0.0, 1.0)  # or adjust if you expect >1
    )
    ax11.legend(loc='upper right', fontsize='small')

    # 7) Save the 2×2 figure
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    # 1) Load run charges (only runs ≤ MAX_RUNNUM)
    run_charges = parse_run_charges(CSV_PATH)

    # 2) Make the 2×2 comparison plot
    make_four_panel_plots(NH3_FILES, C_FILES, H2_FILES, run_charges, FOUR_PANEL_OUTPUT)


if __name__ == "__main__":
    main()