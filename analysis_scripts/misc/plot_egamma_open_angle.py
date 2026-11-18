
import uproot
import numpy as np
import matplotlib.pyplot as plt

files = {
    "CLASDIS Fa22": "/volatile/clas12/thayward/egamma_study/clasdis_fa22_egamma.root",
    "RGC Fa22": "/volatile/clas12/thayward/egamma_study/rgc_fa22_egamma.root",
}

plt.figure(figsize=(10, 6))

bins = np.linspace(0, 180, 181)

for label, path in files.items():
    with uproot.open(path) as f:
        angles = f["PhysicsEvents"]["open_angle"].array(library="np")
    #endwith

    angles = angles[np.isfinite(angles)]
    angles = angles[(angles >= 0) & (angles <= 180)]

    plt.hist(
        angles,
        bins=bins,
        histtype="step",
        linewidth=1.8,
        density=True,
        label=f"{label} (N={len(angles):,})",
    )
#endfor

plt.xlabel(r"Electron-photon opening angle (deg)")
plt.ylabel("Normalized counts")
plt.title(r"CLAS12 Fa22: Electron-photon opening angle")
plt.xlim(0, 180)
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()

output = "/u/home/thayward/egamma_open_angle.png"
plt.savefig(output, dpi=200)
plt.close()

print(f"Saved {output}")
