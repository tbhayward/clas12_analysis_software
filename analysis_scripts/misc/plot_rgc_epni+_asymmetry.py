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
# HARD‐CODED DILUTION DATA (for bottom‐right overlay)
# -----------------------------------------------------------------------------

# Su22 dilution points
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

# Fa22 dilution points
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

# Sp23 dilution points (if/when available)
enpichi2FitsALUsinphi_Sp23 = [
    [0.091314751, 0.060825979, 0.081001615],
    [0.166606555, 0.070263513, 0.012697034],
    [0.251219557, 0.107175756, 0.009614248],
    [0.346388218, 0.141913308, 0.010857952],
    [0.441235552, 0.118073067, 0.014919450],
    [0.535187978, 0.115828125, 0.026600146]
]
x_Sp23      = np.array([row[0] for row in enpichi2FitsALUsinphi_Sp23])
dil_Sp23    = np.array([0.0]*6)       # Placeholder; fill in when ready
dil_err_Sp23= np.array([0.0]*6)

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
    Runs outside these → Eb = 0.
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
    Compute Mandelstam t event‐by‐event for ep → e' n π⁺, assuming exclusivity:
    t = (q - p_π)², where q = p_beam - p_e'.

    Arguments are identical‐length NumPy arrays:
      run_arr, e_p_arr, e_th_arr, e_ph_arr, p_p_arr, p_th_arr, p_ph_arr
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
    Worker for parallel histogramming in Mx².  Returns:
      (label, hist_Mx2_t, counts_base, counts_t)
    where:
      hist_Mx2_t = charge‐normalized Mx² histogram (|t|<1),
      counts_base = number of events surviving run/Mx² ≤2.0 filter,
      counts_t    = number of events surviving also |t|<1.
    Both histograms use 100 bins on [0,1.5].
    """
    filepath, label, run_charges = args

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
    mask_base  = (run_arr <= MAX_RUNNUM) & (mx2_arr <= 2.0)
    run_base   = run_arr[mask_base]
    mx2_base   = mx2_arr[mask_base]
    e_p_base   = e_p_arr[mask_base]
    e_th_base  = e_th_arr[mask_base]
    e_ph_base  = e_ph_arr[mask_base]
    p_p_base   = p_p_arr[mask_base]
    p_th_base  = p_th_arr[mask_base]
    p_ph_base  = p_ph_arr[mask_base]

    # Number of events after base filter
    num_base_events = run_base.size
    print(f"[DEBUG][{label}] after base filter: {num_base_events} events")

    # If no events survive base filter, return zeros
    bins = np.linspace(0.0, 1.5, 101)
    hist_Mx2_t = np.zeros(len(bins)-1, dtype=float)
    if num_base_events == 0:
        return (label, hist_Mx2_t, 0, 0)

    # Compute t_array on base‐filtered events
    t_vals = compute_t_array(
        run_base,
        e_p_base, e_th_base, e_ph_base,
        p_p_base, p_th_base, p_ph_base
    )
    mask_t    = np.abs(t_vals) < 1.0
    run_t     = run_base[mask_t]
    mx2_t     = mx2_base[mask_t]

    num_t_events = run_t.size
    print(f"[DEBUG][{label}] after |t|<1 cut: {num_t_events} events")

    # Build charge‐normalized |t|<1 Mx² histogram
    if num_t_events > 0:
        unique_runs_t = set(np.unique(run_t))
        Q_t = 0.0
        for r in unique_runs_t:
            q = run_charges.get(r, 0.0)
            if q <= 0.0:
                print(f"    [WARNING] (t) run {r} has zero or missing charge.")
            Q_t += q

        if Q_t > 0.0 and mx2_t.size > 0:
            counts_t, _   = np.histogram(mx2_t, bins=bins)
            hist_Mx2_t    = counts_t.astype(float) / Q_t

    return (label, hist_Mx2_t, num_base_events, num_t_events)


def make_four_panel_plot(nh3_files, c_files, h2_files, run_charges, outpath):
    """
    Build a 2×2 figure:
      ┌─────────────────────────┬─────────────────────────┐
      │ (0,0): NH₃ + H₂         │ (0,1): C + H₂           │
      │         (|t|<1, Mx²≤2)  │         (|t|<1, Mx²≤2)  │
      ├─────────────────────────┼─────────────────────────┤
      │ (1,0): (NH₃ – C)/HQ_n   │ (1,1): [ (H_n – s·H_c) / H_n ] (in x_B) │
      │          (|t|<1)        │        (0.75<Mx²<1.05,|t|<1)            │
      └─────────────────────────┴────────────────────────────────────────────┘

    Here:
     • H_n(Mx²): NH₃ charge‐normalized Mx² histogram (|t|<1).
     • H_c(Mx²): C   charge‐normalized Mx² histogram (|t|<1).
     • scale := ∫_{0→0.5} H_n(Mx²)  /  ∫_{0→0.5} H_c(Mx²).
     • Panel (1,1) reuses that same “scale” and plots [H_n(x_B) – scale·H_c(x_B)]/H_n(x_B),
       using events restricted to 0.75<Mx²<1.05, |t|<1.
    """

    # 1) In parallel, build charge‐normalized Mx² histograms (|t|<1) and event counts
    tasks = []
    for fp, lbl in nh3_files + c_files + h2_files:
        tasks.append((fp, lbl, run_charges))

    results = {}  # label → (hist_Mx2_t, num_base, num_t)
    with ProcessPoolExecutor() as executor:
        future_to_label = {executor.submit(process_file, args): args[1] for args in tasks}
        for future in as_completed(future_to_label):
            lbl = future_to_label[future]
            try:
                label, hist_Mx2_t, num_base, num_t = future.result()
                results[label] = (hist_Mx2_t, num_base, num_t)
            except Exception as e:
                print(f"[ERROR] processing {lbl}: {e}")
                results[lbl] = (np.zeros(100, dtype=float), 0, 0)

    # 2) Compute Mx²‐based scale for each period (using bins 0→0.5)
    bins_Mx2 = np.linspace(0.0, 1.5, 101)
    bin_centers_Mx2 = 0.5 * (bins_Mx2[:-1] + bins_Mx2[1:])
    mask_0_5 = (bin_centers_Mx2 >= 0.0) & (bin_centers_Mx2 < 0.5)

    period_names = ["Su22", "Fa22", "Sp23"]
    scale_dict = {}
    for period in period_names:
        h_n = results[f"{period}-NH3"][0]
        h_c = results[f"{period}-C"][0]
        sum_n = h_n[mask_0_5].sum()
        sum_c = h_c[mask_0_5].sum()
        scale = (sum_n / sum_c) if (sum_c > 0.0) else 0.0
        scale_dict[period] = scale
        print(f"[DEBUG] scale({period}) = {scale:.4f}  [∑H_n(0-0.5)={sum_n:.4f}, ∑H_c(0-0.5)={sum_c:.4f}]")

    # 3) Now build x_B histograms (|t|<1, 0.75<Mx²<1.05) for NH₃ and C (per period)
    #    We do this inside a small helper loop (no parallelization needed because NH₃/C are just 3+3 files).
    def build_xB_hist(filepath, run_charges):
        """
        Returns (hist_xB_t, edges_xB) for events satisfying:
          run ≤ MAX_RUNNUM, 0.75<Mx²<1.05, |t|<1.
        Normalized by total charge from unique runs passing that selection.
        """
        tree = uproot.open(filepath)["PhysicsEvents"]
        run_arr  = tree["runnum"].array(library="np").astype(int)
        mx2_arr  = tree["Mx2"].array(library="np")
        e_p_arr  = tree["e_p"].array(library="np")
        e_th_arr = tree["e_theta"].array(library="np")
        e_ph_arr = tree["e_phi"].array(library="np")
        p_p_arr  = tree["p_p"].array(library="np")
        p_th_arr = tree["p_theta"].array(library="np")
        p_ph_arr = tree["p_phi"].array(library="np")
        xB_arr   = tree["x"].array(library="np")

        # 1) Base filter: run ≤ MAX_RUNNUM & 0.75 < Mx² < 1.05
        mask_base_xB = (run_arr <= MAX_RUNNUM) & (mx2_arr > 0.75) & (mx2_arr < 1.05)
        run_base_xB  = run_arr[mask_base_xB]
        mx2_base_xB  = mx2_arr[mask_base_xB]
        e_p_base_xB  = e_p_arr[mask_base_xB]
        e_th_base_xB = e_th_arr[mask_base_xB]
        e_ph_base_xB = e_ph_arr[mask_base_xB]
        p_p_base_xB  = p_p_arr[mask_base_xB]
        p_th_base_xB = p_th_arr[mask_base_xB]
        p_ph_base_xB = p_ph_arr[mask_base_xB]
        xB_base      = xB_arr[mask_base_xB]

        if run_base_xB.size == 0:
            return np.zeros(100, dtype=float), np.linspace(0.0,1.0,101)

        # 2) Compute t and mask |t|<1
        t_vals_xB = compute_t_array(
            run_base_xB,
            e_p_base_xB, e_th_base_xB, e_ph_base_xB,
            p_p_base_xB, p_th_base_xB, p_ph_base_xB
        )
        mask_t_xB = np.abs(t_vals_xB) < 1.0
        run_t_xB  = run_base_xB[mask_t_xB]
        xB_t      = xB_base[mask_t_xB]

        if run_t_xB.size == 0:
            return np.zeros(100, dtype=float), np.linspace(0.0,1.0,101)

        # 3) Sum charge over unique runs in this selection
        uniq_runs_t_xB = set(np.unique(run_t_xB))
        Q_t_xB = sum(run_charges.get(r, 0.0) for r in uniq_runs_t_xB)

        # 4) Build histogram in x_B (we choose 100 bins on [0,1])
        edges_xB = np.linspace(0.0, 1.0, 101)
        counts_xB, _ = np.histogram(xB_t, bins=edges_xB)

        if Q_t_xB > 0.0:
            hist_xB_t = counts_xB.astype(float) / Q_t_xB
        else:
            hist_xB_t = np.zeros(len(edges_xB)-1, dtype=float)

        return hist_xB_t, edges_xB

    # Build x_B histograms for each NH₃/C file
    xB_hist_n = {}
    xB_hist_c = {}
    edges_xB   = None
    for (fp, lbl) in nh3_files:
        hist, edges = build_xB_hist(fp, run_charges)
        xB_hist_n[lbl] = hist
        edges_xB = edges
    for (fp, lbl) in c_files:
        hist, _      = build_xB_hist(fp, run_charges)
        xB_hist_c[lbl] = hist

    # 4) Begin plotting: 2 rows × 2 columns
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharey=False)

    # ──────────────────────────────────────────────────────────────────────────────
    # Panel (0,0): NH₃ + H₂ (|t|<1, Mx²≤2)
    # ──────────────────────────────────────────────────────────────────────────────
    top0_vals = []
    for _, lbl in NH3_FILES:
        h_n = results[lbl][0]
        axes[0, 0].step(bins_Mx2[:-1], h_n, where='post', label=lbl)
        top0_vals.append(h_n.max())

    h2_lbl = H2_FILES[0][1]
    h2_curve = results[h2_lbl][0]
    axes[0, 0].step(
        bins_Mx2[:-1], h2_curve, where='post',
        color='k', linestyle='-', label=h2_lbl
    )
    top0_vals.append(h2_curve.max())

    axes[0, 0].set(
        title="NH₃ + H₂: Normalized Mₓ² (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        ylabel="events / nC",
        xlim=(0.0, 1.5)
    )
    # dynamic y‐limit
    y0 = 1.1 * max(top0_vals) if top0_vals else 0.1
    axes[0, 0].set(ylim=(0.0, y0))
    axes[0, 0].legend(loc='upper right', fontsize='small')

    # ──────────────────────────────────────────────────────────────────────────────
    # Panel (0,1): C + H₂ (|t|<1, Mx²≤2)
    # ──────────────────────────────────────────────────────────────────────────────
    top1_vals = []
    for _, lbl in C_FILES:
        h_c = results[lbl][0]
        axes[0, 1].step(bins_Mx2[:-1], h_c, where='post', label=lbl)
        top1_vals.append(h_c.max())

    axes[0, 1].step(
        bins_Mx2[:-1], h2_curve, where='post',
        color='k', linestyle='-', label=h2_lbl
    )
    top1_vals.append(h2_curve.max())

    axes[0, 1].set(
        title="C + H₂: Normalized Mₓ² (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    y1 = 1.1 * max(top1_vals) if top1_vals else 0.1
    axes[0, 1].set(ylim=(0.0, y1))
    axes[0, 1].legend(loc='upper right', fontsize='small')

    # ──────────────────────────────────────────────────────────────────────────────
    # Panel (1,0): (NH₃ – C) using same Mx² histograms (|t|<1)
    # ──────────────────────────────────────────────────────────────────────────────
    bot0_vals = []
    for period in period_names:
        h_n = results[f"{period}-NH3"][0]
        h_c = results[f"{period}-C"][0]
        scale = scale_dict[period]
        diff = h_n - scale * h_c
        style = ['-', '--', '-.', ':'][period_names.index(period) % 4]
        axes[1, 0].step(
            bins_Mx2[:-1], diff, where='post',
            linestyle=style, label=f"{period} NH₃–C"
        )
        bot0_vals.append(diff.max())

    axes[1, 0].step(
        bins_Mx2[:-1], h2_curve, where='post',
        color='k', linestyle='-', label=h2_lbl
    )
    bot0_vals.append(h2_curve.max())

    axes[1, 0].set(
        title="(NH₃ – C) + H₂ (|t|<1)",
        xlabel=r"$M_x^2$ (GeV$^2$)",
        xlim=(0.0, 1.5)
    )
    y2 = 1.1 * max(bot0_vals) if bot0_vals else 0.1
    axes[1, 0].set(ylim=(0.0, y2))
    axes[1, 0].legend(loc='upper right', fontsize='small')

    # ──────────────────────────────────────────────────────────────────────────────
    # Panel (1,1): [H_n(x_B) – scale·H_c(x_B)] / H_n(x_B) (|t|<1, 0.75<Mx²<1.05)
    # ──────────────────────────────────────────────────────────────────────────────
    bot1_vals = []
    for period in period_names:
        h_n_xB = xB_hist_n[f"{period}-NH3"]
        h_c_xB = xB_hist_c[f"{period}-C"]
        scale  = scale_dict[period]
        # Avoid division by zero: mask bins where H_n(xB)==0
        ratio = np.zeros_like(h_n_xB)
        mask_nonzero = h_n_xB > 0
        ratio[mask_nonzero] = (h_n_xB[mask_nonzero] - scale * h_c_xB[mask_nonzero]) / h_n_xB[mask_nonzero]
        style = ['-', '--', '-.', ':'][period_names.index(period) % 4]
        axes[1, 1].step(
            edges_xB[:-1], ratio, where='post',
            linestyle=style, label=f"{period} (Hₙ–s·H_c)/Hₙ"
        )
        bot1_vals.append(np.nanmax(ratio))

    # Overlay dilution‐factor points in same color
    # Su22:
    axes[1, 1].errorbar(
        x_Su22, dil_Su22, yerr=dil_err_Su22,
        fmt='o', color=axes[1, 1].lines[0].get_color(), label="Su22 dilution"
    )
    # Fa22:
    axes[1, 1].errorbar(
        x_Fa22, dil_Fa22, yerr=dil_err_Fa22,
        fmt='o', color=axes[1, 1].lines[1].get_color(), label="Fa22 dilution"
    )
    # Sp23 (if data exists):
    if np.any(dil_Sp23):
        axes[1, 1].errorbar(
            x_Sp23, dil_Sp23, yerr=dil_err_Sp23,
            fmt='o', color=axes[1, 1].lines[2].get_color(), label="Sp23 dilution"
        )

    axes[1, 1].set(
        title=r"(H$_n$–s·H$_c$)/H$_n$ vs. $x_B$ (|t|<1, 0.75<M$_x^2$<1.05)",
        xlabel=r"$x_B$",
        xlim=(0.0, 1.0)
    )
    y3 = 1.1 * max(bot1_vals + [np.nanmax(dil_Su22), np.nanmax(dil_Fa22), np.nanmax(dil_Sp23)]) if bot1_vals else 1.0
    axes[1, 1].set(ylim=(0.0, y3))
    axes[1, 1].legend(loc='upper right', fontsize='small')

    # 5) Finalize & save
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

    # 2) Create the 2×2 comparison plot
    make_four_panel_plot(NH3_FILES, C_FILES, H2_FILES, run_charges, FOUR_PANEL_OUTPUT)


if __name__ == "__main__":
    main()