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
    [-1.146121838, 0.054859991, 0.012520412],
    [-0.946174495, 0.062743625, 0.010920577],
    [-0.746158937, 0.045566561, 0.009683464],
    [-0.544204731, 0.063586464, 0.008309587],
    [-0.340337870, 0.087657371, 0.006787094],
    [-0.126605312, 0.097653216, 0.005160013],
]

rga_midlow = [
    [-1.145897758, 0.099462884, 0.010818572],
    [-0.945638230, 0.105271301, 0.009566552],
    [-0.744436635, 0.116331259, 0.008286868],
    [-0.542901059, 0.118921909, 0.006979998],
    [-0.342058348, 0.118082449, 0.005922069],
    [-0.133472843, 0.121517711, 0.004417913],
]

rga_midhigh = [
    [-1.145457763, 0.121193538, 0.010745241],
    [-0.945170595, 0.150203716, 0.009437912],
    [-0.744388646, 0.140809516, 0.008196579],
    [-0.544265379, 0.142506153, 0.007123189],
    [-0.343381748, 0.142313217, 0.006180098],
    [-0.141577336, 0.138493603, 0.005046466],
]

rga_high = [
    [-1.146263114, 0.171469353, 0.012038303],
    [-0.945900160, 0.154026517, 0.010856988],
    [-0.745568037, 0.135166243, 0.009665831],
    [-0.545493191, 0.136636233, 0.008591928],
    [-0.344952399, 0.145066824, 0.007505281],
    [-0.150500294, 0.115813850, 0.006880107],
]

# RGC
rgc_low = [
    [-1.145517323, 0.054590352, 0.029159931],
    [-0.945524007, 0.037122901, 0.025406287],
    [-0.745335939, 0.043202871, 0.022268795],
    [-0.543347006, 0.077626797, 0.018945470],
    [-0.340057539, 0.130790439, 0.015112879],
    [-0.126050052, 0.104729304, 0.010064230],
]

rgc_midlow = [
    [-1.145292910, 0.158573764, 0.024681249],
    [-0.945709339, 0.059566901, 0.021898473],
    [-0.743794625, 0.115902149, 0.018933353],
    [-0.543108099, 0.114807548, 0.015625096],
    [-0.341028106, 0.164865490, 0.012740278],
    [-0.134615059, 0.130033768, 0.009010673],
]

rgc_midhigh = [
    [-1.145124545, 0.137028682, 0.024032565],
    [-0.944767425, 0.108568083, 0.021228654],
    [-0.744884316, 0.149548280, 0.018325235],
    [-0.543808893, 0.133456776, 0.015702187],
    [-0.343168094, 0.152778117, 0.013226314],
    [-0.145344286, 0.125388886, 0.010690175],
]

rgc_high = [
    [-1.146819715, 0.143241768, 0.026392794],
    [-0.945992484, 0.154467403, 0.023548312],
    [-0.745673749, 0.106509499, 0.020713445],
    [-0.545265148, 0.120385376, 0.018181117],
    [-0.344533969, 0.151682077, 0.015625239],
    [-0.155687425, 0.093515548, 0.014893490],
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
rga_pull = np.zeros_like(rga_y)
rgc_pull = np.zeros_like(rga_y)

total_chi2 = 0.0

for i in range(len(rga_y)):
    w_rga = 1.0 / (rga_err[i] * rga_err[i])
    w_rgc = 1.0 / (rgc_err[i] * rgc_err[i])

    mean_y[i] = (w_rga * rga_y[i] + w_rgc * rgc_y[i]) / (w_rga + w_rgc)

    rga_pull[i] = (rga_y[i] - mean_y[i]) / rga_err[i]
    rgc_pull[i] = (rgc_y[i] - mean_y[i]) / rgc_err[i]

    total_chi2 += ((rga_y[i] - mean_y[i]) / rga_err[i]) ** 2
    total_chi2 += ((rgc_y[i] - mean_y[i]) / rgc_err[i]) ** 2
#endfor

# 48 total measurements - 24 fitted means
ndf_total = len(rga_y)
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
ax_top.set_title(r"Exclusive $\pi^{+}$ comparison: RGA vs RGC", fontsize=15)

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
ax_bot.set_ylim(-8.0, 8.0)
ax_bot.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("output/exclusive_pip_rga_rgc_comparison.png", dpi=300)
plt.close()

print("Saved plot to: output/exclusive_pip_rga_rgc_comparison.png")