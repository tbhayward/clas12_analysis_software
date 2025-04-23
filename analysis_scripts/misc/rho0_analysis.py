import os
import uproot                              # pip install uproot awkward
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

# === CUT LABELS ===
# no missing-mass cuts
CUT_LABEL_BASE = (
    "$Q^{2}>2, W>2, y<0.75$\n"
    "$z_{\\rho}>0.9, \\Delta\\pi^{-}<0.05$"
)
# with missing-mass cuts
CUT_LABEL_MX = (
    "$Q^{2}>2, W>2, y<0.75$\n"
    "$z_{\\rho}>0.9, |M_{x}^{2}|<0.01$\n"
    "$M_{x\\pi^{+}}^{2}>3.24, M_{x\\pi^{-}}^{2}>1.8225, \\Delta\\pi^{-}<0.05$"
)

# === 1) KINEMATIC‐CUTS FUNCTION ===
def kinematic_cuts(events, use_missing_mass_cuts=False, beam_energy=10.6):
    """
    Return a boolean mask for events passing the base cuts.
    If use_missing_mass_cuts=True, also apply the tighter Mx² cuts.
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
    px_miss = - (px_e  + px_p1 + px_p2)
    py_miss = - (py_e  + py_p1 + py_p2)
    pz_miss = beam_energy - (pz_e + pz_p1 + pz_p2)

    p_miss = np.sqrt(px_miss**2 + py_miss**2 + pz_miss**2)
    p_miss[p_miss == 0] = np.nan  # avoid division by zero

    theta_calc  = np.arccos(pz_miss / p_miss)
    delta_theta = np.abs(p3_theta - theta_calc)

    # base cuts
    mask = (
        (Q2    > 2.0)    &
        (W     > 2.0)    &
        (y     < 0.75)   &
        (z_rho > 0.9)    &
        (delta_theta < 0.05)
    )

    if use_missing_mass_cuts:
        mask &= (
            (np.abs(Mx2)   < 0.01) &
            (Mx2_2  > 3.24)        &
            (Mx2_3  > 1.8225)
        )
    #endif

    return mask
    #endif


# === 2) ORIGINAL MISSING‐MASK PLOTTING FUNCTION ===
def plot_missing_masses(data_dict, mask_dict):
    """
    1×3 panel for Mx2, Mx2_2, Mx2_3 without missing-mass cuts in mask,
    each with its own y-axis scale and x-limits.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    specs = [
        {"key": "Mx2",    "xlabel": r"$M_{x}^{2}\,(\mathrm{GeV}^2)$", 
         "xlim": (-0.03, 0.03)},
        {"key": "Mx2_2",  "xlabel": r"$M_{x\pi^{+}}^{2}\,(\mathrm{GeV}^2)$", 
         "xlim": (1, 6)},
        {"key": "Mx2_3",  "xlabel": r"$M_{x\pi^{-}}^{2}\,(\mathrm{GeV}^2)$", 
         "xlim": (1, 6)},
    ]

    for ax, spec in zip(axes, specs):
        # collect masked arrays
        arrs = [ev[spec["key"]][mask_dict[label]]
                for label, ev in data_dict.items()]
        labels = list(data_dict.keys())

        # per-panel y_max
        max_count = 0
        hists = []
        for arr in arrs:
            counts, edges = np.histogram(arr, bins=100, range=spec["xlim"])
            hists.append((counts, edges))
            max_count = max(max_count, counts.max())
        y_max = max_count * 1.2 if max_count > 0 else 1.0

        # plot
        for (counts, edges), label in zip(hists, labels):
            centers = 0.5 * (edges[:-1] + edges[1:])
            errors  = np.sqrt(counts)
            ax.errorbar(centers, counts, yerr=errors, fmt='o', label=label)

        ax.set_xlabel(spec["xlabel"])
        ax.set_ylabel("Counts")
        ax.set_xlim(*spec["xlim"])
        ax.set_ylim(0, y_max)

        # annotate base cuts
        ax.text(
            0.05, 0.95, CUT_LABEL_BASE,
            transform=ax.transAxes, va='top', fontsize=8
        )
        ax.legend(loc='upper right', fontsize=8)
    #endfor

    plt.tight_layout()
    os.makedirs("output/rho0", exist_ok=True)
    plt.savefig("output/rho0/missing_masses.pdf")
    plt.close()


# === 3) BEAM-SPIN ASYMMETRY CALCULATION ===
def calculate_bsa(data_dict):
    """
    Compute BSA vs cos(theta) in 16 bins from -1 to 1.
    In each cosθ bin, histogram phi23 (12 bins 0→2π), fit A*sinϕ,
    then scale by <DepA>/(<beam_pol>*<DepW>).
    Returns dict[label] = (cos_centers, bsa_vals, bsa_errs).
    Also returns a combined average.
    """
    cos_bins = np.linspace(-1, 1, 17)
    phi_bins = np.linspace(0, 2*np.pi, 13)

    results = {}
    for label, ev in data_dict.items():
        # apply full cuts including missing-mass
        mask = kinematic_cuts(ev, use_missing_mass_cuts=True,
                               beam_energy=BEAM_ENERGIES[label])
        theta = ev["theta"][mask]
        cos_t = np.cos(theta)
        phi23   = ev["phi23"][mask]
        DepA    = ev["DepA"][mask]
        DepW    = ev["DepW"][mask]
        pol     = ev["beam_pol"][mask]

        cos_centers = 0.5 * (cos_bins[:-1] + cos_bins[1:])
        bsa_vals = np.zeros_like(cos_centers)
        bsa_errs = np.zeros_like(cos_centers)

        for i in range(len(cos_centers)):
            idx = (cos_t >= cos_bins[i]) & (cos_t < cos_bins[i+1])
            if idx.sum() < 3:
                bsa_vals[i] = np.nan
                bsa_errs[i] = np.nan
                continue
            # histogram phi23
            counts, edges = np.histogram(phi23[idx], bins=phi_bins)
            centers = 0.5 * (edges[:-1] + edges[1:])
            x = np.sin(centers)
            y = counts

            # least-squares fit A*sinϕ
            A = np.sum(x*y) / np.sum(x*x)
            res = y - A*x
            var_res = np.sum(res**2) / (len(y)-1) if len(y)>1 else 0
            sigma_A = np.sqrt(var_res/np.sum(x*x)) if np.sum(x*x)>0 else np.nan

            # scale by mean DepA/(pol*DepW)
            scale = np.mean(DepA[idx]) / (np.mean(pol[idx]) * np.mean(DepW[idx]))
            bsa_vals[i] = A * scale
            bsa_errs[i] = sigma_A * scale

        results[label] = (cos_centers, bsa_vals, bsa_errs)
    #endfor

    # combined average
    cos_centers = cos_centers
    vs = np.array([results[label][1] for label in results])
    es = np.array([results[label][2] for label in results])
    w = np.where(es>0, 1/es**2, 0)
    num = np.nansum(w * vs, axis=0)
    den = np.nansum(w, axis=0)
    comb = np.where(den>0, num/den, np.nan)
    comb_err = np.where(den>0, np.sqrt(1/den), np.nan)
    results["Combined"] = (cos_centers, comb, comb_err)

    return results


# === 4) PLOTTING BSA ===
def plot_bsa(bsa_dict):
    """
    1×2 panel: left = individual datasets, right = combined.
    Annotate full cuts.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

    # left: each dataset
    for label in DATASETS.keys():
        if label in bsa_dict:
            x, y, yerr = bsa_dict[label]
            axes[0].errorbar(x, y, yerr=yerr, fmt='o', label=label)
    axes[0].text(0.05, 0.95, CUT_LABEL_MX,
                 transform=axes[0].transAxes, va='top', fontsize=8)
    axes[0].set_xlabel(r"$\cos\theta$")
    axes[0].set_ylabel(r"$F_{LU}/F_{UU}$")
    axes[0].set_xlim(-1, 1)
    axes[0].set_ylim(-1, 0.2)
    axes[0].legend(loc='upper right', fontsize=8)

    # right: combined
    x, y, yerr = bsa_dict["Combined"]
    axes[1].errorbar(x, y, yerr=yerr, fmt='o', label="Combined")
    axes[1].text(0.05, 0.95, CUT_LABEL_MX,
                 transform=axes[1].transAxes, va='top', fontsize=8)
    axes[1].set_xlabel(r"$\cos\theta$")
    axes[1].set_xlim(-1, 1)
    axes[1].set_ylim(-1, 0.2)
    axes[1].legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    os.makedirs("output/rho0", exist_ok=True)
    plt.savefig("output/rho0/BSA.pdf")
    plt.close()


# === 5) MAIN DRIVER ===
def main():
    data_dict = {}
    mask_dict = {}

    # load and mask (base cuts) for missing-mass plot
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
            # for missing-mass plot: no Mx² cuts
            mask_dict[label] = kinematic_cuts(arrs,
                                              use_missing_mass_cuts=False,
                                              beam_energy=BEAM_ENERGIES[label])
    #endfor

    # make the missing-mass panel
    plot_missing_masses(data_dict, mask_dict)

    # now compute & plot BSA (with missing-mass cuts)
    bsa = calculate_bsa(data_dict)
    plot_bsa(bsa)


if __name__ == "__main__":
    main()
    #endif