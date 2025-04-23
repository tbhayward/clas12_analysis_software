import os
import uproot                              # pip install uproot awkward
import numpy as np
import matplotlib.pyplot as plt

# === CONFIGURATION ===
QUICK_RUN   = False                        # set False to loop over full tree
MAX_EVENTS  = 100_000                     # only used when QUICK_RUN is True

DATASETS = {
    # "RGA Fa18 Inb": (
    #     "/work/clas12/thayward/CLAS12_SIDIS/processed_data/"
    #     "pass2/data/eppi+pi-X/rga_fa18_inb_eppi+pi-X.root"
    # ),
    # "RGA Fa18 Out": (
    #     "/work/clas12/thayward/CLAS12_SIDIS/processed_data/"
    #     "pass2/data/eppi+pi-X/rga_fa18_out_eppi+pi-X.root"
    # ),
    # "RGA Sp19 Inb": (
    #     "/work/clas12/thayward/CLAS12_SIDIS/processed_data/"
    #     "pass2/data/eppi+pi-X/rga_sp19_inb_eppi+pi-X.root"
    # ),
    "RGA Fa18 Inb": (
        "/scratch/thayward/"
        "rga_fa18_inb_eppi+pi-X.root"
    ),
    "RGA Fa18 Out": (
        "/scratch/thayward/"
        "rga_fa18_out_eppi+pi-X.root"
    ),
    "RGA Sp19 Inb": (
        "/scratch/thayward/"
        "rga_sp19_inb_eppi+pi-X.root"
    ),
}

# beam energies by dataset (GeV)
BEAM_ENERGIES = {
    "RGA Fa18 Inb": 10.6,
    "RGA Fa18 Out": 10.6,
    "RGA Sp19 Inb": 10.2,
}

TREE_NAME = "PhysicsEvents"

# === CUT LABELS ===
CUT_LABEL_BASE = (
    "$Q^{2}>2, W>2, y<0.75$\n"
    "$z_{\\rho}>0.9, \\Delta\\pi^{-}<0.05$"
)
CUT_LABEL_MX = (
    "$Q^{2}>2, W>2, y<0.75$\n"
    "$z_{\\rho}>0.9, |M_{x}^{2}|<0.01$\n"
    "$M_{x\\pi^{+}}^{2}>3.24, M_{x\\pi^{-}}^{2}>1.8225, \\Delta\\pi^{-}<0.05$"
)

# === 1) KINEMATIC‐CUTS FUNCTION ===
def kinematic_cuts(ev, use_missing_mass_cuts=False, beam_energy=10.6):
    Q2, W, y, z_rho = ev["Q2"], ev["W"], ev["y"], ev["z23"]
    Mx2, Mx2_2, Mx2_3 = ev["Mx2"], ev["Mx2_2"], ev["Mx2_3"]
    p3_theta = ev["p3_theta"]
    e_p, e_theta, e_phi = ev["e_p"], ev["e_theta"], ev["e_phi"]
    p1_p, p1_theta, p1_phi = ev["p1_p"], ev["p1_theta"], ev["p1_phi"]
    p2_p, p2_theta, p2_phi = ev["p2_p"], ev["p2_theta"], ev["p2_phi"]

    def to_cart(p, th, ph):
        return (p * np.sin(th) * np.cos(ph),
                p * np.sin(th) * np.sin(ph),
                p * np.cos(th))

    px_e, py_e, pz_e    = to_cart(e_p,  e_theta,  e_phi)
    px_p1, py_p1, pz_p1 = to_cart(p1_p, p1_theta, p1_phi)
    px_p2, py_p2, pz_p2 = to_cart(p2_p, p2_theta, p2_phi)

    px_m = -(px_e  + px_p1 + px_p2)
    py_m = -(py_e  + py_p1 + py_p2)
    pz_m = beam_energy - (pz_e + pz_p1 + pz_p2)

    p_m = np.sqrt(px_m**2 + py_m**2 + pz_m**2)
    p_m[p_m == 0] = np.nan
    theta_calc  = np.arccos(pz_m / p_m)
    delta_th    = np.abs(p3_theta - theta_calc)

    # mask = (
    #     (Q2    > 2.0) &
    #     (W     > 2.0) &
    #     (y     < 0.75) &
    #     (z_rho > 0.9) &
    #     (delta_th < 0.05)
    # )
    # if use_missing_mass_cuts:
    #     mask &= (
    #         (np.abs(Mx2) < 0.01) &
    #         (Mx2_2     > 3.24) &
    #         (Mx2_3     > 1.8225)
    #     )
    mask = (
        (Q2    > 0.0) &
        (W     > 0.0) &
        (y     < 1) &
        (z_rho > 0) &
        (delta_th < 50)
    )
    if use_missing_mass_cuts:
        mask &= (
            (np.abs(Mx2) <100) &
            (Mx2_2     > 0) &
            (Mx2_3     > 0)
    )
    return mask
    #endif


# === 2) MISSING‐MASS PLOTTING ===
def plot_missing_masses(data, masks):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    specs = [
        ("Mx2",   r"$M_{x}^{2}\,(\mathrm{GeV}^2)$",    (-0.03, 0.03)),
        ("Mx2_2", r"$M_{x\pi^{+}}^{2}\,(\mathrm{GeV}^2)$",(1, 6)),
        ("Mx2_3", r"$M_{x\pi^{-}}^{2}\,(\mathrm{GeV}^2)$",(1, 6)),
    ]
    for ax, (key, xlabel, xlim) in zip(axes, specs):
        ymax = 0
        for label, ev in data.items():
            arr = ev[key][masks[label]]
            counts, edges = np.histogram(arr, bins=100, range=xlim)
            centers = 0.5*(edges[:-1]+edges[1:])
            ax.errorbar(centers, counts, yerr=np.sqrt(counts),
                        fmt='o', label=label)
            ymax = max(ymax, counts.max())
        ymax = ymax*1.2 if ymax>0 else 1.0

        ax.set(xlabel=xlabel, ylabel="Counts", xlim=xlim, ylim=(0, ymax))
        ax.text(0.05, 0.95, CUT_LABEL_BASE, transform=ax.transAxes,
                va='top', fontsize=8)
        ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    os.makedirs("output/rho0", exist_ok=True)
    plt.savefig("output/rho0/missing_masses.pdf")
    plt.close()


# === 3) BEAM‐SPIN ASYMMETRY CALCULATION ===
def calculate_bsa(data):
    cos_bins = np.linspace(-1, 1, 17)
    phi_bins = np.linspace(0, 2*np.pi, 13)
    results = {}

    for label, ev in data.items():
        mask = kinematic_cuts(ev, use_missing_mass_cuts=True,
                               beam_energy=BEAM_ENERGIES[label])
        cos_t   = np.cos(ev["theta"][mask])
        phi23   = ev["phi23"][mask]
        helic   = ev["helicity"][mask]
        DepA    = ev["DepA"][mask]
        DepW    = ev["DepW"][mask]
        pol     = ev["beam_pol"][mask]

        centers = 0.5*(cos_bins[:-1] + cos_bins[1:])
        vals = np.zeros_like(centers)
        errs = np.zeros_like(centers)

        for i in range(len(centers)):
            sel = (cos_t >= cos_bins[i]) & (cos_t < cos_bins[i+1])
            if sel.sum() < 3:
                vals[i] = np.nan
                errs[i] = np.nan
                continue

            Np, edges = np.histogram(phi23[sel & (helic > 0)], bins=phi_bins)
            Nm, _     = np.histogram(phi23[sel & (helic < 0)], bins=phi_bins)
            tot = Np + Nm
            with np.errstate(divide='ignore', invalid='ignore'):
                asym = (Np - Nm) / tot

            ph_centers = 0.5*(edges[:-1] + edges[1:])
            x = np.sin(ph_centers)
            y = asym

            A = np.nansum(x*y) / np.nansum(x*x)
            res = y - A*x
            var = np.nansum(res**2)/(len(y)-1) if len(y)>1 else 0
            sigmaA = np.sqrt(var/np.nansum(x*x)) if np.nansum(x*x)>0 else np.nan

            scale = np.mean(DepA[sel]) / (np.mean(pol[sel])*np.mean(DepW[sel]))
            vals[i] = A * scale
            errs[i] = sigmaA * scale

        results[label] = (centers, vals, errs)

    # combined
    all_vals = np.array([results[l][1] for l in results])
    all_errs = np.array([results[l][2] for l in results])
    w = np.where(all_errs>0, 1/all_errs**2, 0)
    num = np.nansum(all_vals*w, axis=0)
    den = np.nansum(w, axis=0)
    comb = np.where(den>0, num/den, np.nan)
    comb_err = np.where(den>0, np.sqrt(1/den), np.nan)
    centers = next(iter(results.values()))[0]
    results["Combined"] = (centers, comb, comb_err)

    return results


# === 4) PLOTTING BSA ===
def plot_bsa(bsa):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

    for label in DATASETS:
        if label in bsa:
            x, y, ye = bsa[label]
            axes[0].errorbar(x, y, yerr=ye, fmt='o', label=label)
    axes[0].text(0.05, 0.95, CUT_LABEL_MX, transform=axes[0].transAxes,
                 va='top', fontsize=8)
    axes[0].set(xlabel=r"$\cos\theta$", ylabel=r"$F_{LU}/F_{UU}$",
                xlim=(-1, 1), ylim=(-1, 1))
    axes[0].legend(loc='upper right', fontsize=8)

    x, y, ye = bsa["Combined"]
    axes[1].errorbar(x, y, yerr=ye, fmt='o', label="Combined")
    axes[1].text(0.05, 0.95, CUT_LABEL_MX, transform=axes[1].transAxes,
                 va='top', fontsize=8)
    axes[1].set(xlabel=r"$\cos\theta$", xlim=(-1, 1), ylim=(-1, 1))
    axes[1].legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    os.makedirs("output/rho0", exist_ok=True)
    plt.savefig("output/rho0/BSA.pdf")
    plt.close()


# === 5) θ AND cosθ DISTRIBUTIONS ===
def plot_theta_distributions(data, masks):
    """
    1×2 panel: left = θ distribution, right = cosθ distribution,
    using base kinematic cuts.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    bins_theta = np.linspace(0, np.pi, 50)
    bins_cos   = np.linspace(-1, 1, 50)

    for label, ev in data.items():
        mask = masks[label]
        theta = ev["theta"][mask]
        cos_t = np.cos(theta)
        axes[0].hist(theta, bins=bins_theta, histtype='step', label=label)
        axes[1].hist(cos_t, bins=bins_cos, histtype='step', label=label)
    #endfor

    axes[0].set(xlabel=r"$\theta$ (rad)", ylabel="Counts")
    axes[1].set(xlabel=r"$\cos\theta$", ylabel="Counts")
    axes[0].legend(loc='upper right', fontsize=8)
    axes[1].legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    os.makedirs("output/rho0", exist_ok=True)
    plt.savefig("output/rho0/theta_distributions.pdf")
    plt.close()


# === 6) MAIN DRIVER ===
def main():
    data = {}
    masks = {}
    branches = [
        "Q2","W","y","z23",
        "Mx2","Mx2_2","Mx2_3",
        "p3_theta",
        "e_p","e_theta","e_phi",
        "p1_p","p1_theta","p1_phi",
        "p2_p","p2_theta","p2_phi",
        "theta","phi23","helicity","DepA","DepW","beam_pol"
    ]

    for label, path in DATASETS.items():
        with uproot.open(path) as f:
            tree = f[TREE_NAME]
            stop = MAX_EVENTS if QUICK_RUN else None
            arrs = tree.arrays(branches, entry_stop=stop, library="np")
            data[label] = arrs
            masks[label] = kinematic_cuts(
                arrs, use_missing_mass_cuts=False,
                beam_energy=BEAM_ENERGIES[label]
            )
    #endfor

    plot_missing_masses(data, masks)
    bsa = calculate_bsa(data)
    plot_bsa(bsa)
    plot_theta_distributions(data, masks)


if __name__ == "__main__":
    main()
    #endif