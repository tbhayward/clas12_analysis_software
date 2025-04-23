import os
import uproot                              # pip install uproot
import numpy as np
import matplotlib.pyplot as plt

# === CONFIGURATION ===
QUICK_RUN   = True                        # set False to loop over the full tree
MAX_EVENTS  = 100_000                     # only used when QUICK_RUN is True

DATASETS = {
    "RGA Fa18 Inb": (
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/"
        "pass2/data/eppi+pi-X/rga_fa18_inb_eppi+pi-X.root"
    ),
    # "RGA Fa18 Out": (
    #     "/work/clas12/thayward/CLAS12_SIDIS/processed_data/"
    #     "pass2/data/eppi+pi-X/rga_fa18_out_eppi+pi-X.root"
    # ),
    "RGA Sp19 Inb": (
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/"
        "pass2/data/eppi+pi-X/rga_sp19_inb_eppi+pi-X.root"
    ),
}

# beam energies by dataset (GeV)
BEAM_ENERGIES = {
    "RGA Fa18 Inb": 10.6,
    # "RGA Fa18 Out": 10.6,
    "RGA Sp19 Inb": 10.2,
}

TREE_NAME = "PhysicsEvents"

# This LaTeX‐style string will be drawn in the top‐left of each panel,
# split across three lines and in smaller text.
CUT_LABEL = (
    r"$Q^{2}>2,\ W>2,\ y<0.75$\\" + "\n"
    r"$z_{\rho}>0.9,\ |M_{x}^{2}|<0.3$\\" + "\n"
    r"$M_{x\pi^{+}}^{2}>3.24,\ M_{x\pi^{-}}^{2}>1.83,\ \Delta\pi^{-}<0.05$"
)

# === 1) KINEMATIC‐CUTS FUNCTION ===
def kinematic_cuts(events, use_missing_mass_cuts=False, beam_energy=10.6):
    """
    Given an 'events' dict of NumPy arrays for one dataset,
    return a boolean mask where all selected cuts pass.
    If use_missing_mass_cuts=True, also apply the three Mx² cuts.
    beam_energy (GeV) is used to compute the missing π⁻ direction.
    """
    # extract branches
    Q2        = events["Q2"]
    W         = events["W"]
    y         = events["y"]
    z_rho     = events["z23"]
    Mx2       = events["Mx2"]
    Mx2_2     = events["Mx2_2"]
    Mx2_3     = events["Mx2_3"]
    p3_theta  = events["p3_theta"]

    e_p       = events["e_p"]
    e_theta   = events["e_theta"]
    e_phi     = events["e_phi"]
    p1_p      = events["p1_p"]
    p1_theta  = events["p1_theta"]
    p1_phi    = events["p1_phi"]
    p2_p      = events["p2_p"]
    p2_theta  = events["p2_theta"]
    p2_phi    = events["p2_phi"]

    # helper: convert (p,θ,φ) to cartesian
    def to_cart(p, theta, phi):
        px = p * np.sin(theta) * np.cos(phi)
        py = p * np.sin(theta) * np.sin(phi)
        pz = p * np.cos(theta)
        return px, py, pz

    px_e,  py_e,  pz_e   = to_cart(e_p,  e_theta,  e_phi)
    px_p1, py_p1, pz_p1 = to_cart(p1_p, p1_theta, p1_phi)
    px_p2, py_p2, pz_p2 = to_cart(p2_p, p2_theta, p2_phi)

    # missing momentum = beam + target - (e'+p'+π+)
    px_miss = - (px_e +   px_p1 + px_p2)
    py_miss = - (py_e +   py_p1 + py_p2)
    pz_miss = beam_energy - (pz_e + pz_p1 + pz_p2)

    p_miss  = np.sqrt(px_miss**2 + py_miss**2 + pz_miss**2)
    p_miss[p_miss == 0] = np.nan  # avoid div by zero

    theta_calc   = np.arccos(pz_miss / p_miss)
    delta_theta  = np.abs(p3_theta - theta_calc)

    # base mask (all cuts except missing‐mass)
    mask = (
        (Q2    > 2.0)    &
        (W     > 2.0)    &
        (y     < 0.75)   &
        (z_rho > 0.9)    &
        (delta_theta < 0.05)
    )

    if use_missing_mass_cuts:
        mask &= (
            (np.abs(Mx2)   < 0.3) &
            (Mx2_2 > 3.24)       &
            (Mx2_3 > 1.8225)
            # & (Mx2_extra1 > ... )
            # & (Mx2_extra2 > ... )
        )
    #endif

    return mask
    #endif


# === 2) PLOTTING FUNCTION ===
def plot_missing_masses(data_dict, mask_dict):
    """
    Build a 1×3 panel for Mx2, Mx2_2, Mx2_3 distributions
    from each dataset in data_dict, using their masks in mask_dict.
    """
    Mxs, Mx2s, Mx3s, labels = [], [], [], []
    for label, ev in data_dict.items():
        m = mask_dict[label]
        Mxs .append(ev["Mx2"][m])
        Mx2s.append(ev["Mx2_2"][m])
        Mx3s.append(ev["Mx2_3"][m])
        labels.append(label)
    #endfor

    # compute common y‐axis max
    all_max = []
    for arr_list in (Mxs, Mx2s, Mx3s):
        for arr in arr_list:
            counts, _ = np.histogram(arr, bins=100)
            all_max.append(counts.max())
    y_max = max(all_max) * 1.2
    if y_max <= 0:
        y_max = 1.0

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    x_labels = [
        r"$M_{x}^{2}\,\mathrm{(GeV^{2})}$",
        r"$M_{x\pi^{+}}^{2}\,\mathrm{(GeV^{2})}$",
        r"$M_{x\pi^{-}}^{2}\,\mathrm{(GeV^{2})}$",
    ]

    for ax, arr_list, xlabel in zip(axes, (Mxs, Mx2s, Mx3s), x_labels):
        for arr, label in zip(arr_list, labels):
            counts, edges = np.histogram(arr, bins=100, range=(-1, 4))
            centers = 0.5*(edges[:-1] + edges[1:])
            errors  = np.sqrt(counts)
            ax.errorbar(centers, counts, yerr=errors, fmt='o', label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Counts")
        ax.set_xlim(-1, 4)
        ax.set_ylim(0, y_max)
    #endfor

    axes[0].text(
        0.05, 0.95, CUT_LABEL,
        transform=axes[0].transAxes,
        va='top', fontsize=8
    )
    axes[2].legend(loc='upper right')

    plt.tight_layout()
    os.makedirs("output/rho0", exist_ok=True)
    plt.savefig("output/rho0/missing_masses.pdf")
    plt.close()


# === 3) MAIN DRIVER ===
def main():
    data_dict = {}
    mask_dict = {}

    for label, path in DATASETS.items():
        with uproot.open(path) as f:
            tree = f[TREE_NAME]
            stop = MAX_EVENTS if QUICK_RUN else None

            arrs = tree.arrays(
                [
                    "Q2","W","y","z23",
                    "Mx2","Mx2_2","Mx2_3",
                    "p3_theta",
                    "e_p","e_theta","e_phi",
                    "p1_p","p1_theta","p1_phi",
                    "p2_p","p2_theta","p2_phi"
                ],
                entry_stop=stop,
                library="np"
            )

            data_dict[label] = arrs
            beam_E = BEAM_ENERGIES[label]
            # for initial plot: exclude missing‐mass cuts
            mask_dict[label] = kinematic_cuts(
                arrs,
                use_missing_mass_cuts=False,
                beam_energy=beam_E
            )
    #endfor

    plot_missing_masses(data_dict, mask_dict)


if __name__ == "__main__":
    main()
    #endif