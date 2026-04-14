import os
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

# ------------------------------------------------
# input data
# each triplet is:
# [tprime_value, bsa_value, uncertainty]
# ------------------------------------------------

# RGA
rga_low = [
    [-1.156284951, 0.046264860, 0.000776900],
    [-0.956016087, 0.063803215, 0.000967777],
    [-0.755285835, 0.072591575, 0.001186778],
    [-0.554449590, 0.075764803, 0.001423349],
    [-0.353922308, 0.084719413, 0.001673759],
    [-0.150698758, 0.101250471, 0.001880106],
]

rga_midlow = [
    [-1.156284951, 0.046264860, 0.000776900],
    [-0.956016087, 0.063803215, 0.000967777],
    [-0.755285835, 0.072591575, 0.001186778],
    [-0.554449590, 0.075764803, 0.001423349],
    [-0.353922308, 0.084719413, 0.001673759],
    [-0.150698758, 0.101250471, 0.001880106],
]

rga_midhigh = [
    [-1.156284951, 0.046264860, 0.000776900],
    [-0.956016087, 0.063803215, 0.000967777],
    [-0.755285835, 0.072591575, 0.001186778],
    [-0.554449590, 0.075764803, 0.001423349],
    [-0.353922308, 0.084719413, 0.001673759],
    [-0.150698758, 0.101250471, 0.001880106],
]

rga_high = [
    [-1.156284951, 0.046264860, 0.000776900],
    [-0.956016087, 0.063803215, 0.000967777],
    [-0.755285835, 0.072591575, 0.001186778],
    [-0.554449590, 0.075764803, 0.001423349],
    [-0.353922308, 0.084719413, 0.001673759],
    [-0.150698758, 0.101250471, 0.001880106],
]

# RGC
rgc_low = [
    [-1.156625652, 0.046189715, 0.001273375],
    [-0.956165924, 0.066534536, 0.001594316],
    [-0.755566163, 0.080916683, 0.001960521],
    [-0.555411688, 0.087594526, 0.002378914],
    [-0.355860968, 0.108280348, 0.002871685],
    [-0.158412692, 0.097091866, 0.003586773],
]

rgc_midlow = [
    [-1.156284951, 0.046264860, 0.000776900],
    [-0.956016087, 0.063803215, 0.000967777],
    [-0.755285835, 0.072591575, 0.001186778],
    [-0.554449590, 0.075764803, 0.001423349],
    [-0.353922308, 0.084719413, 0.001673759],
    [-0.150698758, 0.101250471, 0.001880106],
]

rgc_midhigh = [
    [-1.156284951, 0.046264860, 0.000776900],
    [-0.956016087, 0.063803215, 0.000967777],
    [-0.755285835, 0.072591575, 0.001186778],
    [-0.554449590, 0.075764803, 0.001423349],
    [-0.353922308, 0.084719413, 0.001673759],
    [-0.150698758, 0.101250471, 0.001880106],
]

rgc_high = [
    [-1.156284951, 0.046264860, 0.000776900],
    [-0.956016087, 0.063803215, 0.000967777],
    [-0.755285835, 0.072591575, 0.001186778],
    [-0.554449590, 0.075764803, 0.001423349],
    [-0.353922308, 0.084719413, 0.001673759],
    [-0.150698758, 0.101250471, 0.001880106],
]

# ------------------------------------------------
# helpers
# ------------------------------------------------

def flatten_reversed(groups):
    values = []
    errors = []

    for group in groups:
        reversed_group = list(reversed(group))
        for point in reversed_group:
            values.append(point[1])
            errors.append(point[2])
        #endfor
    #endfor

    return np.array(values, dtype=float), np.array(errors, dtype=float)
#enddef

# ------------------------------------------------
# flatten into 24 bins
# ------------------------------------------------

rga_groups = [rga_low, rga_midlow, rga_midhigh, rga_high]
rgc_groups = [rgc_low, rgc_midlow, rgc_midhigh, rgc_high]

rga_y, rga_err = flatten_reversed(rga_groups)
rgc_y, rgc_err = flatten_reversed(rgc_groups)

x = np.arange(1, len(rga_y) + 1, dtype=float)

if len(rga_y) != 24 or len(rgc_y) != 24:
    raise RuntimeError("Expected 24 points in each dataset.")
#endif

# ------------------------------------------------
# weighted mean bin-by-bin, pulls, chi2
# ------------------------------------------------

mean_y = np.zeros_like(rga_y)
mean_err = np.zeros_like(rga_y)

rga_pull = np.zeros_like(rga_y)
rgc_pull = np.zeros_like(rga_y)

total_chi2 = 0.0
ndf = len(rga_y)

for i in range(len(rga_y)):
    w_rga = 1.0 / (rga_err[i] * rga_err[i])
    w_rgc = 1.0 / (rgc_err[i] * rgc_err[i])

    mean_y[i] = (w_rga * rga_y[i] + w_rgc * rgc_y[i]) / (w_rga + w_rgc)
    mean_err[i] = math.sqrt(1.0 / (w_rga + w_rgc))

    rga_pull[i] = (rga_y[i] - mean_y[i]) / rga_err[i]
    rgc_pull[i] = (rgc_y[i] - mean_y[i]) / rgc_err[i]

    total_chi2 += ((rga_y[i] - mean_y[i]) / rga_err[i]) ** 2
    total_chi2 += ((rgc_y[i] - mean_y[i]) / rgc_err[i]) ** 2
#endfor

# There are 48 measurements and 24 fitted means
ndf_total = 2 * len(rga_y) - len(rga_y)
chi2_per_ndf = total_chi2 / ndf_total
p_value = chi2.sf(total_chi2, ndf_total)

print("")
print("============================================================")
print("RGA vs RGC consistency test about common weighted mean")
print("============================================================")
print(f"Number of bins           : {len(rga_y)}")
print(f"Total measurements       : {2 * len(rga_y)}")
print(f"Degrees of freedom       : {ndf_total}")
print(f"Total chi2               : {total_chi2:.6f}")
print(f"chi2/ndf                 : {chi2_per_ndf:.6f}")
print(f"p-value                  : {p_value:.6e}")
print("============================================================")
print("")

# ------------------------------------------------
# plot
# ------------------------------------------------

os.makedirs("output", exist_ok=True)

fig, (ax_top, ax_bot) = plt.subplots(
    2,
    1,
    figsize=(14, 9),
    sharex=True,
    gridspec_kw={"height_ratios": [2.2, 1.0]}
)

# top panel
ax_top.errorbar(
    x - 0.08,
    rga_y,
    yerr=rga_err,
    fmt="o",
    markersize=5,
    capsize=3,
    label="RGA"
)

ax_top.errorbar(
    x + 0.08,
    rgc_y,
    yerr=rgc_err,
    fmt="s",
    markersize=5,
    capsize=3,
    label="RGC"
)

ax_top.set_ylabel(r"$F_{LU}^{\sin\phi}/F_{UU}$", fontsize=14)
ax_top.legend(fontsize=12)
ax_top.grid(True, alpha=0.3)

# bottom panel
ax_bot.axhline(0.0, linestyle="--", linewidth=1.0)
ax_bot.axhline(1.0, linestyle=":", linewidth=1.0)
ax_bot.axhline(-1.0, linestyle=":", linewidth=1.0)
ax_bot.axhline(2.0, linestyle=":", linewidth=1.0)
ax_bot.axhline(-2.0, linestyle=":", linewidth=1.0)

ax_bot.plot(
    x - 0.08,
    rga_pull,
    "o",
    markersize=5,
    label="RGA pull"
)

ax_bot.plot(
    x + 0.08,
    rgc_pull,
    "s",
    markersize=5,
    label="RGC pull"
)

ax_bot.set_xlabel("bin number", fontsize=14)
ax_bot.set_ylabel("pull", fontsize=14)
ax_bot.set_xlim(0.5, 24.5)
ax_bot.grid(True, alpha=0.3)

title = (
    r"Exclusive $\pi^{+}$ comparison: RGA vs RGC"
    + "\n"
    + rf"$\chi^2/\mathrm{{ndf}} = {chi2_per_ndf:.3f}$, "
    + rf"$p = {p_value:.3e}$"
)
ax_top.set_title(title, fontsize=15)

plt.tight_layout()
plt.savefig("output/exclusive_pip_rga_rgc_comparison.png", dpi=300)
plt.close()

print("Saved plot to: output/exclusive_pip_rga_rgc_comparison.png")