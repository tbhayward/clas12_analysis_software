#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

# Toggle for quick debugging (process only first 5 runs per period)
QUICK_RUN = True

# Maximum run number to include
MAX_RUNNUM = 17768

# Path to the CSV of run charges
CSV_PATH = "/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

# Output directory and filenames
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

m_e    = 0.000511   # electron mass [GeV]
m_pi   = 0.13957    # charged pion mass [GeV]
m_p    = 0.938272   # proton mass [GeV]
m_n    = 0.939565   # neutron mass [GeV]  (not strictly needed below, but included for clarity)

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

    with open(csv_path, 'r') as f:
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


def load_tree(filepath):
    """
    Open a ROOT file and return the 'PhysicsEvents' TTree.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Could not find file {filepath}")
    tree = uproot.open(filepath)["PhysicsEvents"]
    print(f"[DEBUG] Loaded TTree 'PhysicsEvents' from {filepath} with {tree.num_entries} entries.")
    return tree


def get_beam_energy(runnums):
    """
    Given a NumPy array of run numbers, return a NumPy array of beam energies Eb (in GeV) event-by-event.
    - 16042 ≤ run ≤ 17065 → Eb = 10.5473 GeV
    - 17067 ≤ run ≤ 17724 → Eb = 10.5563 GeV
    - 17725 ≤ run ≤ 17811 → Eb = 10.5593 GeV
    If runnum falls outside these ranges, Eb = 0 (will be masked out later).
    """
    Eb = np.zeros_like(runnums, dtype=float)
    mask1 = (runnums >= 16042) & (runnums <= 17065)
    mask2 = (runnums >= 17067) & (runnums <= 17724)
    mask3 = (runnums >= 17725) & (runnums <= 17811)
    Eb[mask1] = 10.5473
    Eb[mask2] = 10.5563
    Eb[mask3] = 10.5593
    return Eb


def compute_t_array(run_arr, e_p_arr, e_theta_arr, e_phi_arr, p_p_arr, p_theta_arr, p_phi_arr):
    """
    Compute the Mandelstam variable t (for the neutron) on an event-by-event basis,
    assuming ep → e' n π⁺ is exclusive.  We do this via:

      q := p_beam - p_e'   (virtual photon 4-vector)
      p_pi := 4-momentum of the detected π⁺
      t = (q - p_pi)²

    where (q - p_pi)² = (E_q - E_pi)² - |⃗q - ⃗p_pi|².

    - run_arr:       array of run numbers (so we know Eb for each event).
    - e_p_arr:       array of scattered-electron momentum magnitudes [GeV].
    - e_theta_arr:   polar angle of e' [radians].
    - e_phi_arr:     azimuthal angle of e' [radians].
    - p_p_arr:       array of π⁺ momentum magnitudes [GeV].
    - p_theta_arr:   polar angle of π⁺ [radians].
    - p_phi_arr:     azimuthal angle of π⁺ [radians].

    Returns:
      t_array: 1D NumPy array of t (in GeV²), same length as inputs.
    """
    # 1) Determine Eb for each event:
    Eb_arr = get_beam_energy(run_arr)

    # 2) Build 4-vector for scattered electron e':
    #    E_e' = sqrt(|⃗p_e'|² + m_e²)
    E_e_arr = np.sqrt(e_p_arr**2 + m_e**2)
    #    p_e'_x = e_p sinθ_e cosφ_e,  p_e'_y = e_p sinθ_e sinφ_e,  p_e'_z = e_p cosθ_e
    sin_e = np.sin(e_theta_arr)
    cos_e = np.cos(e_theta_arr)
    ex = e_p_arr * sin_e * np.cos(e_phi_arr)
    ey = e_p_arr * sin_e * np.sin(e_phi_arr)
    ez = e_p_arr * cos_e

    # 3) Build 4-vector for π⁺:
    #    E_π = sqrt(|⃗p_π|² + m_pi²)
    E_pi_arr = np.sqrt(p_p_arr**2 + m_pi**2)
    sin_p = np.sin(p_theta_arr)
    cos_p = np.cos(p_theta_arr)
    px = p_p_arr * sin_p * np.cos(p_phi_arr)
    py = p_p_arr * sin_p * np.sin(p_phi_arr)
    pz = p_p_arr * cos_p

    # 4) Beam 4-vector p_beam = (Eb, 0, 0, +Eb)  (electron traveling along +z)
    # 5) Virtual photon q = p_beam - p_e'
    #    => E_q = Eb - E_e'
    #       q_x = 0 - ex = -ex
    #       q_y = 0 - ey = -ey
    #       q_z = Eb - ez
    E_q  = Eb_arr - E_e_arr
    qx   = -ex
    qy   = -ey
    qz   = Eb_arr - ez

    # 6) (q - p_π) 4-vector = (E_q - E_π,  q⃗ - p⃗_π)
    E_diff = E_q - E_pi_arr
    dx     = qx - px
    dy     = qy - py
    dz     = qz - pz

    # 7) Compute t = (E_diff)² - (dx² + dy² + dz²)
    t_array = (E_diff**2) - (dx**2 + dy**2 + dz**2)
    return t_array


def compute_charge_and_mask(tree, run_charges, quick=False):
    """
    For a given TTree, return (Q_tot, masked_Mx2_array) where:
      - Q_tot is the sum of charges for unique runs ≤ MAX_RUNNUM
        (only first 5 if quick=True),
      - masked_Mx2_array are the Mx² values for those runs
        with Mx² ≤ 2.0 *and* |t| < 1.0.

    In other words, we apply three filters:
      1) runnum ≤ MAX_RUNNUM
      2) Mx² ≤ 2.0
      3) |t| < 1.0

    We extract all relevant arrays (runnum, Mx², e_p, e_theta, e_phi, p_p, p_theta, p_phi)
    and compute t for each event, then mask accordingly.
    """
    # Read runnum, Mx², and kinematic branches:
    run_arr    = tree["runnum"].array(library="np").astype(int)
    mx2_arr    = tree["Mx2"].array(library="np")
    e_p_arr    = tree["e_p"].array(library="np")
    e_th_arr   = tree["e_theta"].array(library="np")
    e_ph_arr   = tree["e_phi"].array(library="np")
    p_pi_arr   = tree["p_p"].array(library="np")
    p_th_arr   = tree["p_theta"].array(library="np")
    p_ph_arr   = tree["p_phi"].array(library="np")

    # 1) Filter out runs > MAX_RUNNUM AND Mx² > 2.0
    base_mask = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
    run_filt  = run_arr[base_mask]
    mx2_filt  = mx2_arr[base_mask]
    e_p_filt  = e_p_arr[base_mask]
    e_th_filt = e_th_arr[base_mask]
    e_ph_filt = e_ph_arr[base_mask]
    p_p_filt  = p_pi_arr[base_mask]
    p_th_filt = p_th_arr[base_mask]
    p_ph_filt = p_ph_arr[base_mask]

    if run_filt.size == 0:
        return 0.0, np.array([], dtype=float)

    # 2) Compute t-array for the filtered events
    t_vals = compute_t_array(
        run_filt,
        e_p_filt, e_th_filt, e_ph_filt,
        p_p_filt, p_th_filt, p_ph_filt
    )

    # 3) Now apply |t| < 1.0
    t_masked = np.abs(t_vals) < 1.0

    # Combine masks:
    run_final = run_filt[t_masked]
    mx2_final = mx2_filt[t_masked]

    if run_final.size == 0:
        return 0.0, np.array([], dtype=float)

    # 4) If QUICK_RUN, select only the first 5 unique runnums encountered
    if quick:
        seen = []
        for r in run_final:
            if r not in seen:
                seen.append(r)
            if len(seen) >= 5:
                break
        selected_runs = set(seen)
    else:
        selected_runs = set(np.unique(run_final))

    # 5) Sum total charge over those selected runs
    Q_tot = 0.0
    for r in selected_runs:
        q = run_charges.get(r, 0.0)
        if q <= 0.0:
            print(f"    [WARNING] run {r} has zero or missing charge.")
        Q_tot += q

    # 6) Build final mask (to keep only events from selected_runs)
    final_mask = np.isin(run_final, list(selected_runs))
    return Q_tot, mx2_final[final_mask]


def make_normalized_Mx2_plots(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 1×3 figure:
      - Panel 1: NH₃ (Su22, Fa22, Sp23) + H₂
      - Panel 2:  C   (Su22, Fa22, Sp23) + H₂
      - Panel 3: Differences [NH₃ – C] for each period (normalized to match 0–0.5), plus H₂

    All histograms apply these cuts:
      • runnum ≤ MAX_RUNNUM
      • Mx² ≤ 2.0
      • |t| < 1.0

    (If QUICK_RUN=True, only the first 5 unique runs per file are used.)
    """
    bins = np.linspace(0.0, 1.5, 101)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    # -------------------------------------------------------------
    # Panel 1: NH₃ + H₂
    # -------------------------------------------------------------
    for filepath, label in nh3_files:
        tree = load_tree(filepath)
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0.0 or mx2_vals.size == 0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        norm_counts = counts.astype(float) / Q
        axes[0].step(bins[:-1], norm_counts, where='post', label=label)

    for filepath, label in h2_files:
        tree = load_tree(filepath)
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0.0 or mx2_vals.size == 0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        norm_counts = counts.astype(float) / Q
        axes[0].step(
            bins[:-1], norm_counts, where='post',
            color='k', linestyle='-', label=label
        )

    axes[0].set(
        title="NH₃ + H₂: Normalized Mx² (|t|<1)",
        xlabel=r"$M_x^2$ [GeV$^2$]",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    axes[0].legend(loc='upper right', fontsize='small')

    # -------------------------------------------------------------
    # Panel 2: C + H₂
    # -------------------------------------------------------------
    for filepath, label in c_files:
        tree = load_tree(filepath)
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0.0 or mx2_vals.size == 0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        norm_counts = counts.astype(float) / Q
        axes[1].step(bins[:-1], norm_counts, where='post', label=label)

    for filepath, label in h2_files:
        tree = load_tree(filepath)
        Q, mx2_vals = compute_charge_and_mask(tree, run_charges, quick=QUICK_RUN)
        if Q <= 0.0 or mx2_vals.size == 0:
            continue
        counts, _ = np.histogram(mx2_vals, bins=bins)
        norm_counts = counts.astype(float) / Q
        axes[1].step(
            bins[:-1], norm_counts, where='post',
            color='k', linestyle='-', label=label
        )

    axes[1].set(
        title="C + H₂: Normalized Mx² (|t|<1)",
        xlabel=r"$M_x^2$ [GeV$^2$]",
        xlim=(0.0, 1.5)
    )
    axes[1].legend(loc='upper right', fontsize='small')

    # -------------------------------------------------------------
    # Panel 3: (NH₃ – C) differences, plus H₂
    # -------------------------------------------------------------
    period_names = ["Su22", "Fa22", "Sp23"]
    nh3_hist_dict = {}
    c_hist_dict   = {}

    # 1) Build period‐level histograms for NH₃ and C
    for period_label in period_names:
        # Find the file with that prefix
        nh3_fp = next(fp for fp, lbl in nh3_files if lbl.startswith(period_label))
        c_fp   = next(fp for fp, lbl in c_files   if lbl.startswith(period_label))

        # NH₃:
        tree_n = load_tree(nh3_fp)
        Q_n, mx2_n = compute_charge_and_mask(tree_n, run_charges, quick=QUICK_RUN)
        if Q_n > 0.0 and mx2_n.size > 0:
            counts_n, _ = np.histogram(mx2_n, bins=bins)
            nh3_hist_dict[period_label] = counts_n.astype(float) / Q_n
        else:
            nh3_hist_dict[period_label] = np.zeros(len(bins)-1)

        # C:
        tree_c = load_tree(c_fp)
        Q_c, mx2_c = compute_charge_and_mask(tree_c, run_charges, quick=QUICK_RUN)
        if Q_c > 0.0 and mx2_c.size > 0:
            counts_c, _ = np.histogram(mx2_c, bins=bins)
            c_hist_dict[period_label] = counts_c.astype(float) / Q_c
        else:
            c_hist_dict[period_label] = np.zeros(len(bins)-1)

    # 2) Build H₂ histogram once
    h2_label = h2_files[0][1]
    tree_h2  = load_tree(h2_files[0][0])
    Q_h2, mx2_h2 = compute_charge_and_mask(tree_h2, run_charges, quick=QUICK_RUN)
    if Q_h2 > 0.0 and mx2_h2.size > 0:
        counts_h2, _ = np.histogram(mx2_h2, bins=bins)
        h2_norm = counts_h2.astype(float) / Q_h2
    else:
        h2_norm = np.zeros(len(bins)-1)

    # 3) Compute (NH₃ – scaled C) for each period, matching 0–0.5
    diff_dict = {}
    for period_label in period_names:
        nh3_vals = nh3_hist_dict[period_label]
        c_vals   = c_hist_dict[period_label]

        # Only bins with center in [0, 0.5)
        mask_0_5 = (bin_centers >= 0.0) & (bin_centers < 0.5)
        nh3_sum = nh3_vals[mask_0_5].sum()
        c_sum   = c_vals[mask_0_5].sum()

        if c_sum > 0.0:
            scale = nh3_sum / c_sum
        else:
            scale = 0.0
        print(f"[DEBUG] Period {period_label}: NH₃ sum(0–0.5)={nh3_sum:.3f}, C sum(0–0.5)={c_sum:.3f}, scale={scale:.3f}")

        diff = nh3_vals - (c_vals * scale)
        diff_dict[period_label] = diff

    # 4) Plot the differences and H₂ on panel 3
    linestyles = ['-', '--', '-.', ':']
    for idx, period_label in enumerate(period_names):
        diff = diff_dict[period_label]
        style = linestyles[idx % len(linestyles)]
        axes[2].step(
            bins[:-1], diff, where='post',
            linestyle=style, label=f"{period_label} NH₃–C"
        )

    # Overlay H₂
    axes[2].step(
        bins[:-1], h2_norm, where='post',
        color='k', linestyle='-', label=h2_label
    )

    axes[2].set(
        title="(NH₃ – C) per period, and H₂ (|t|<1)",
        xlabel=r"$M_x^2$ [GeV$^2$]",
        xlim=(0.0, 1.5),
        ylim=(0.0, 0.03)   # lock y-axis at 0.03
    )
    axes[2].legend(loc='upper right', fontsize='small')

    # 5) Save the three-panel figure
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    # Load run charges, skipping runnum > MAX_RUNNUM
    run_charges = parse_run_charges(CSV_PATH)

    # Create the three‐panel comparison plot (all with |t| < 1)
    make_normalized_Mx2_plots(NH3_FILES, C_FILES, H2_FILES, run_charges, THREE_PANEL_OUTPUT)

if __name__ == "__main__":
    main()