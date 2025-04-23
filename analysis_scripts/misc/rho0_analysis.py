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

# Multiline cut label (no "\\" artifacts)
CUT_LABEL = (
    "$Q^{2}>2, W>2, y<0.75$\n"
    "$z_{\\rho}>0.9, |M_{x}^{2}|<0.3$\n"
    "$M_{x\\pi^{+}}^{2}>3.24, M_{x\\pi^{-}}^{2}>1.83, \\Delta\\pi^{-}<0.05$"
)

# === 1) KINEMATIC‐CUTS FUNCTION ===
def kinematic_cuts(events, use_missing_mass_cuts=False, beam_energy=10.6):
    """
    Return a boolean mask for events passing the base cuts.
    If use_missing_mass_cuts=True, also apply the Mx² cuts.
    """
    Q2       = events["Q2"]
    W        = events["W"]
    y        = events["y"]
    z_rho    = events["z23"]
    Mx2      = events["Mx2"]
    Mx2_2    = events["Mx2_2"]
    Mx2_3    = events["Mx2_3"]
    p3_theta = events["p3_theta"]

    e_p      = events["e_p"]
    e_theta  = events["e_theta"]
    e_phi    = events["e_phi"]
    p1_p     = events["p1_p"]
    p1_theta = events["p1_theta"]
    p1_phi   = events["p1_phi"]
    p2_p     = events["p2_p"]
    p2_theta = events["p2_theta"]
    p2_phi   = events["p2_phi"]

    def to_cart(p, theta, phi):
        px = p * np.sin(theta) * np.cos(phi)
        py = p * np.sin(theta) * np.sin(phi)
        pz = p * np.cos(theta)
        return px, py, pz

    px_e, py_e, pz_e     = to_cart(e_p,  e_theta,  e_phi)
    px_p1, py_p1, pz_p1  = to_cart(p1_p, p1_theta, p1_phi)
    px_p2, py_p2, pz_p2  = to_cart(p2_p, p2_theta, p2_phi)

    px_miss = - (px_e  + px_p1 + px_p2)
    py_miss = - (py_e  + py_p1 + py_p2)
    pz_miss = beam_energy - (pz_e + pz_p1 + pz_p2)

    p_miss = np.sqrt(px_miss**2 + py_miss**2 + pz_miss**2)
    p_miss[p_miss == 0] = np.nan

    theta_calc  = np.arccos(pz_miss / p_miss)
    delta_theta = np.abs(p3_theta - theta_calc)

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
            (Mx2_2  > 3.24)       &
            (Mx2_3  > 1.8225)
        )
    #endif

    return mask
    #endif


# === 2) PLOTTING FUNCTION ===
def plot_missing_masses(data_dict, mask_dict):
    """
    Build a 1×3 panel for Mx2, Mx2_2, Mx2_3 distributions,
    each with its own y-axis scale, x-limits, and dual legend.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    specs = [
        {"key": "Mx2",    "xlabel": r"$M_{x}^{2}\,(\mathrm{GeV}^2)$", 
         "xlim": (-0.05, 0.05)},
        {"key": "Mx2_2",  "xlabel": r"$M_{x\pi^{+}}^{2}\,(\mathrm{GeV}^2)$", 
         "xlim": (1, 6)},
        {"key": "Mx2_3",  "xlabel": r"$M_{x\pi^{-}}^{2}\,(\mathrm{GeV}^2)$", 
         "xlim": (1, 6)},
    ]

    for ax, spec in zip(axes, specs):
        # collect arrays for this quantity
        arrs = [ev[spec["key"]][mask_dict[label]] 
                for label, ev in data_dict.items()]
        labels = list(data_dict.keys())

        # compute per-panel y_max
        max_count = 0
        hist_data = []
        for arr in arrs:
            counts, edges = np.histogram(arr, bins=100, range=spec["xlim"])
            hist_data.append((counts, edges))
            max_count = max(max_count, counts.max())
        y_max = max_count * 1.2 if max_count > 0 else 1.0

        # plot each dataset
        for (counts, edges), label in zip(hist_data, labels):
            centers = 0.5 * (edges[:-1] + edges[1:])
            errors  = np.sqrt(counts)
            ax.errorbar(centers, counts, yerr=errors, fmt='o', label=label)

        ax.set_xlabel(spec["xlabel"])
        ax.set_ylabel("Counts")
        ax.set_xlim(*spec["xlim"])
        ax.set_ylim(0, y_max)

        # text-annotation for cuts
        ax.text(
            0.05, 0.95, CUT_LABEL,
            transform=ax.transAxes, va='top', fontsize=8
        )
        # legend for datasets
        ax.legend(loc='upper right', fontsize=8)
    #endfor

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

            branches = [
                "Q2","W","y","z23",
                "Mx2","Mx2_2","Mx2_3",
                "p3_theta",
                "e_p","e_theta","e_phi",
                "p1_p","p1_theta","p1_phi",
                "p2_p","p2_theta","p2_phi"
            ]
            arrs = tree.arrays(branches, entry_stop=stop, library="np")

            data_dict[label] = arrs
            beam_E = BEAM_ENERGIES[label]
            mask_dict[label] = kinematic_cuts(
                arrs, use_missing_mass_cuts=False, beam_energy=beam_E
            )
    #endfor

    plot_missing_masses(data_dict, mask_dict)


if __name__ == "__main__":
    main()
    #endif