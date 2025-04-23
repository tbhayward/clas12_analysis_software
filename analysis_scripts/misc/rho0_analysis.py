import os
import uproot                              # pip install uproot awkward
import numpy as np
import matplotlib.pyplot as plt

# === CONFIGURATION ===
QUICK_RUN   = True                        # set False to loop over full tree
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

BEAM_ENERGIES = {
    "RGA Fa18 Inb": 10.6,
    # "RGA Fa18 Out": 10.6,
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

# === 1) KINEMATIC‐CUTS ===
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

    px_e, py_e, pz_e   = to_cart(e_p, e_theta, e_phi)
    px_p1, py_p1, pz_p1 = to_cart(p1_p, p1_theta, p1_phi)
    px_p2, py_p2, pz_p2 = to_cart(p2_p, p2_theta, p2_phi)

    px_m = - (px_e + px_p1 + px_p2)
    py_m = - (py_e + py_p1 + py_p2)
    pz_m = beam_energy - (pz_e + pz_p1 + pz_p2)

    p_m = np.sqrt(px_m**2 + py_m**2 + pz_m**2)
    p_m[p_m==0] = np.nan
    theta_calc = np.arccos(pz_m / p_m)
    delta_th = np.abs(p3_theta - theta_calc)

    mask = (
        (Q2    > 2.0) &
        (W     > 2.0) &
        (y     < 0.75) &
        (z_rho > 0.9) &
        (delta_th < 0.05)
    )
    if use_missing_mass_cuts:
        mask &= (
            (np.abs(Mx2) < 0.01) &
            (Mx2_2     > 3.24) &
            (Mx2_3     > 1.8225)
        )
    return mask


# === 2) MISSING‐MASS PLOTS ===
def plot_missing_masses(data, masks):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    specs = [
        ("Mx2",   r"$M_{x}^{2}\,(\mathrm{GeV}^2)$",    (-0.03,0.03)),
        ("Mx2_2", r"$M_{x\pi^{+}}^{2}\,(\mathrm{GeV}^2)$", (1,6)),
        ("Mx2_3", r"$M_{x\pi^{-}}^{2}\,(\mathrm{GeV}^2)$", (1,6)),
    ]
    for ax, (key, xlabel, xlim) in zip(axes, specs):
        hist_data = []
        ymax = 0
        for label, ev in data.items():
            arr = ev[key][masks[label]]
            counts, edges = np.histogram(arr, bins=100, range=xlim)
            hist_data.append((counts, edges, label))
            ymax = max(ymax, counts.max())
        ymax = ymax*1.2 if ymax>0 else 1.0

        for counts, edges, label in hist_data:
            centers = 0.5*(edges[:-1]+edges[1:])
            ax.errorbar(centers, counts, yerr=np.sqrt(counts),
                        fmt='o', label=label)
        ax.set(xlabel=xlabel, ylabel="Counts", xlim=xlim, ylim=(0,ymax))
        ax.text(0.05,0.95, CUT_LABEL_BASE, transform=ax.transAxes,
                va='top', fontsize=8)
        ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    os.makedirs("output/rho0", exist_ok=True)
    plt.savefig("output/rho0/missing_masses.pdf")
    plt.close()


# === 3) BSA CALCULATION ===
def calculate_bsa(data):
    cos_bins = np.linspace(-1,1,17)
    phi_bins = np.linspace(0,2*np.pi,13)
    results = {}

    for label, ev in data.items():
        mask = kinematic_cuts(ev, use_missing_mass_cuts=True,
                               beam_energy=BEAM_ENERGIES[label])
        cos_t = np.cos(ev["theta"][mask])
        phi23 = ev["phi23"][mask]
        DepA  = ev["DepA"][mask]
        DepW  = ev["DepW"][mask]
        pol   = ev["beam_pol"][mask]

        centers = 0.5*(cos_bins[:-1]+cos_bins[1:])
        vals = np.full_like(centers, np.nan)
        errs = np.full_like(centers, np.nan)

        for i in range(len(centers)):
            sel = (cos_t>=cos_bins[i])&(cos_t<cos_bins[i+1])
            if sel.sum()<3: continue
            counts, edges = np.histogram(phi23[sel], bins=phi_bins)
            ph_centers = 0.5*(edges[:-1]+edges[1:])
            x = np.sin(ph_centers)
            y = counts
            A = np.sum(x*y)/np.sum(x*x)
            res = y - A*x
            var = np.sum(res**2)/(len(y)-1) if len(y)>1 else 0
            sigmaA = np.sqrt(var/np.sum(x*x)) if np.sum(x*x)>0 else np.nan
            scale = np.mean(DepA[sel])/(np.mean(pol[sel])*np.mean(DepW[sel]))
            vals[i] = A*scale
            errs[i] = sigmaA*scale

        results[label] = (centers, vals, errs)

    # combined
    all_vals = np.array([results[l][1] for l in results])
    all_errs = np.array([results[l][2] for l in results])
    w = np.where(all_errs>0, 1/all_errs**2, 0)
    num = np.nansum(all_vals*w, axis=0)
    den = np.nansum(w, axis=0)
    comb = np.where(den>0, num/den, np.nan)
    comb_err = np.where(den>0, np.sqrt(1/den), np.nan)
    results["Combined"] = (centers, comb, comb_err)

    return results


# === 4) BSA PLOT ===
def plot_bsa(bsa):
    fig, axes = plt.subplots(1, 2, figsize=(12,5), sharex=True, sharey=True)
    for label in DATASETS:
        if label in bsa:
            x,y,ye = bsa[label]
            axes[0].errorbar(x,y,yerr=ye,fmt='o',label=label)
    axes[0].text(0.05,0.95,CUT_LABEL_MX,transform=axes[0].transAxes,
                 va='top',fontsize=8)
    axes[0].set(xlabel=r"$\cos\theta$",ylabel=r"$F_{LU}/F_{UU}$",
                xlim=(-1,1),ylim=(-1,1))
    axes[0].legend(loc='upper right',fontsize=8)

    x,y,ye = bsa["Combined"]
    axes[1].errorbar(x,y,yerr=ye,fmt='o',label="Combined")
    axes[1].text(0.05,0.95,CUT_LABEL_MX,transform=axes[1].transAxes,
                 va='top',fontsize=8)
    axes[1].set(xlabel=r"$\cos\theta$",xlim=(-1,1),ylim=(-1,1))
    axes[1].legend(loc='upper right',fontsize=8)

    plt.tight_layout()
    plt.savefig("output/rho0/BSA.pdf")
    plt.close()


# === 5) PHI23 FITS DIAGNOSTIC ===
def plot_phi_fits(data):
    cos_bins = np.linspace(-1,1,17)
    phi_bins = np.linspace(0,2*np.pi,13)

    for label, ev in data.items():
        mask = kinematic_cuts(ev, use_missing_mass_cuts=True,
                               beam_energy=BEAM_ENERGIES[label])
        cos_t = np.cos(ev["theta"][mask])
        phi23 = ev["phi23"][mask]

        fig, axes = plt.subplots(4,4, figsize=(12,12), sharex=True)
        axes = axes.flatten()

        for i in range(16):
            ax = axes[i]
            sel = (cos_t>=cos_bins[i])&(cos_t<cos_bins[i+1])
            counts, edges = np.histogram(phi23[sel], bins=phi_bins)
            centers = 0.5*(edges[:-1]+edges[1:])
            ax.errorbar(centers, counts, yerr=np.sqrt(counts), fmt='o')

            # fit
            x = np.sin(centers)
            A = np.sum(x*counts)/np.sum(x*x) if np.sum(x*x)>0 else 0
            phi_line = np.linspace(0,2*np.pi,200)
            ax.plot(phi_line, A*np.sin(phi_line), '-r')

            ax.set_title(f"{centers[0]:.2f}→{centers[-1]:.2f}", fontsize=8)
            ax.set_xlim(0,2*np.pi)

        plt.suptitle(f"Phi23 fits in cosθ bins — {label}")
        plt.tight_layout(rect=[0,0,1,0.96])
        os.makedirs("output/rho0", exist_ok=True)
        fname = label.replace(" ","_")
        plt.savefig(f"output/rho0/phi23_fits_{fname}.pdf")
        plt.close()


# === 6) MAIN DRIVER ===
def main():
    data = {}
    masks = {}
    branches = [
        "Q2","W","y","z23","Mx2","Mx2_2","Mx2_3","p3_theta",
        "e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi",
        "p2_p","p2_theta","p2_phi","theta","phi23","DepA","DepW","beam_pol"
    ]

    for label,path in DATASETS.items():
        with uproot.open(path) as f:
            tree = f[TREE_NAME]
            stop = MAX_EVENTS if QUICK_RUN else None
            arrs = tree.arrays(branches, entry_stop=stop, library="np")
            data[label] = arrs
            masks[label] = kinematic_cuts(arrs,
                                          use_missing_mass_cuts=False,
                                          beam_energy=BEAM_ENERGIES[label])

    plot_missing_masses(data, masks)

    bsa = calculate_bsa(data)
    plot_bsa(bsa)

    plot_phi_fits(data)


if __name__ == "__main__":
    main()