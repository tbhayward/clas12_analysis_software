
import uproot
import numpy as np
import matplotlib.pyplot as plt

files = {
    "CLASDIS Fa22": "/volatile/clas12/thayward/egamma_study/clasdis_fa22_egamma.root",
    "RGC Fa22": "/volatile/clas12/thayward/egamma_study/rgc_fa22_egamma.root",
}

colors = {
    "CLASDIS Fa22": "tab:blue",
    "RGC Fa22": "tab:orange",
}

fig, axes = plt.subplots(2, 1, figsize=(10, 11))

angle_bins = np.linspace(0, 60, 121)
energy_bins = np.linspace(0, 10.6, 107)

for label, path in files.items():
    with uproot.open(path) as f:
        data = f["PhysicsEvents"].arrays(
            ["open_angle", "p_p"], library="np"
        )
    #endwith

    angles = data["open_angle"]
    energies = data["p_p"]

    valid = np.isfinite(angles) & np.isfinite(energies)
    angles = angles[valid]
    energies = energies[valid]

    # Top: opening angle
    angle_mask = (angles >= 0) & (angles <= 60)

    axes[0].hist(
        angles[angle_mask],
        bins=angle_bins,
        histtype="step",
        linewidth=1.8,
        density=True,
        color=colors[label],
        label=f"{label} (N={np.sum(angle_mask):,})",
    )

    # Bottom: photon energy, split at 8 degrees
    for mask, style, region in [
        (angles < 8, "-", r"$\theta_{e\gamma}<8^\circ$"),
        (angles >= 8, "--", r"$\theta_{e\gamma}\geq8^\circ$"),
    ]:
        selection = mask & (angles >= 0) & (angles <= 60)

        if np.any(selection):
            axes[1].hist(
                energies[selection],
                bins=energy_bins,
                histtype="step",
                linewidth=1.8,
                linestyle=style,
                density=True,
                color=colors[label],
                label=f"{label}, {region} (N={np.sum(selection):,})",
            )
        #endif
    #endfor
#endfor

axes[0].set(
    xlabel=r"Electron-photon opening angle (deg)",
    ylabel="Normalized counts",
    title="CLAS12 Fa22: Electron-photon opening angle",
    xlim=(0, 60),
)

axes[1].set(
    xlabel=r"Photon energy (GeV)",
    ylabel="Normalized counts",
    title=r"Photon energy: $\theta_{e\gamma}<8^\circ$ vs. $\geq8^\circ$",
    xlim=(0, 10.6),
)

for ax in axes:
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)
#endfor

plt.tight_layout()

output = "/u/home/thayward/egamma_open_angle.png"
plt.savefig(output, dpi=200)
plt.close()

print(f"Saved {output}")