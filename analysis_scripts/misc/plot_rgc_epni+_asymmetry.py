#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Hard-coded data for all three run periods
#
# Updated Su22, Fa22, and Sp23 data
# -----------------------------------------------------------------------------

# RGC Su22 data (latest fit)
enpichi2FitsALUoffset_Su22 = [
    [0.112750665, 0.071221646, 0.061441005],
    [0.152898421, 0.036786798, 0.012051071],
    [0.198551006, 0.000793036, 0.006367625],
    [0.246623128, 0.017979137, 0.004732863],
    [0.295752215, 0.009408965, 0.004140308],
    [0.344579347, 0.008200952, 0.003995146],
    [0.393290682, 0.007012573, 0.004580434],
    [0.442534108, -0.001417234, 0.005850148],
    [0.491379131, -0.007305088, 0.008338240],
    [0.539671498, -0.007619537, 0.013973362]
]

enpichi2FitsALUsinphi_Su22 = [
    [0.112750665, 0.189086286, 0.118753268],
    [0.152898421, 0.049316026, 0.025837862],
    [0.198551006, 0.094604441, 0.015627219],
    [0.246623128, 0.104457008, 0.012962977],
    [0.295752215, 0.153539915, 0.012586243],
    [0.344579347, 0.142789208, 0.012987987],
    [0.393290682, 0.132304895, 0.014816696],
    [0.442534108, 0.117936342, 0.018100353],
    [0.491379131, 0.167282898, 0.024326647],
    [0.539671498, 0.066769409, 0.038007047]
]

enpichi2FitsAULoffset_Su22 = [
    [0.112750665, 0.017473793, 0.156847694],
    [0.152898421, 0.120767452, 0.030030147],
    [0.198551006, 0.109960519, 0.016957938],
    [0.246623128, 0.100510446, 0.013160185],
    [0.295752215, 0.116720519, 0.013339179],
    [0.344579347, 0.111663937, 0.012972525],
    [0.393290682, 0.097260784, 0.013205055],
    [0.442534108, 0.086136173, 0.015073156],
    [0.491379131, 0.113258956, 0.021011087],
    [0.539671498, 0.121980015, 0.033205902]
]

enpichi2FitsAULsinphi_Su22 = [
    [0.112750665, 0.057332495, 0.130376342],
    [0.152898421, -0.053944105, 0.025308655],
    [0.198551006, -0.039696464, 0.013935427],
    [0.246623128, -0.002814987, 0.010212639],
    [0.295752215, 0.020767347, 0.010204993],
    [0.344579347, 0.030087562, 0.009218074],
    [0.393290682, 0.038931895, 0.010081289],
    [0.442534108, 0.039367852, 0.011867932],
    [0.491379131, 0.055998322, 0.016733904],
    [0.539671498, -0.054528635, 0.027323437]
]

enpichi2FitsAULsin2phi_Su22 = [
    [0.112750665, 0.608918665, 0.376580749],
    [0.152898421, -0.050126078, 0.055950087],
    [0.198551006, -0.002710262, 0.028651610],
    [0.246623128, -0.070282268, 0.022486169],
    [0.295752215, -0.102497822, 0.021434336],
    [0.344579347, -0.102692969, 0.020534720],
    [0.393290682, -0.063108419, 0.021517975],
    [0.442534108, -0.041744960, 0.025523878],
    [0.491379131, -0.050454687, 0.035044107],
    [0.539671498, -0.056709529, 0.054550886]
]

enpichi2FitsALL_Su22 = [
    [0.112750665, -0.727862941, 0.278571679],
    [0.152898421, 0.180118699, 0.052272028],
    [0.198551006, 0.242509048, 0.034657113],
    [0.246623128, 0.318612165, 0.033024802],
    [0.295752215, 0.415067343, 0.038175775],
    [0.344579347, 0.439420451, 0.040484052],
    [0.393290682, 0.538872825, 0.047652324],
    [0.442534108, 0.566609033, 0.053884796],
    [0.491379131, 0.574523329, 0.063127848],
    [0.539671498, 0.619362865, 0.083419778]
]

enpichi2FitsALLcosphi_Su22 = [
    [0.112750665, -0.981175564, 0.501236636],
    [0.152898421, 0.260927358, 0.086720379],
    [0.198551006, 0.216957144, 0.052581937],
    [0.246623128, 0.050988379, 0.049723778],
    [0.295752215, 0.044974519, 0.056032583],
    [0.344579347, -0.022750170, 0.062674287],
    [0.393290682, -0.015125220, 0.071340553],
    [0.442534108, 0.026212836, 0.081782273],
    [0.491379131, -0.151103475, 0.096922172],
    [0.539671498, -0.117320331, 0.131435359]
]

# RGC Fa22 data (latest fit)
enpichi2FitsALUoffset_Fa22 = [
    [0.113905344, -0.003016748, 0.014835380],
    [0.152729530, 0.000370491, 0.008628726],
    [0.198321988, 0.001246482, 0.004701314],
    [0.246523461, -0.000045929, 0.003578731],
    [0.295697442, 0.000117352, 0.003177488],
    [0.344509236, -0.000617232, 0.003105638],
    [0.393376160, -0.002700509, 0.003544395],
    [0.442374288, -0.002261768, 0.004585737],
    [0.491249232, 0.003331318, 0.006628686],
    [0.539365767, -0.000580970, 0.011312685]
]

enpichi2FitsALUsinphi_Fa22 = [
    [0.113905344, -3.439723070, 0.064073166],
    [0.152729530, 0.076281013, 0.018464695],
    [0.198321988, 0.101065977, 0.011471768],
    [0.246523461, 0.125160080, 0.009819273],
    [0.295697442, 0.138006646, 0.009617928],
    [0.344509236, 0.139206399, 0.010114579],
    [0.393376160, 0.140790859, 0.011435606],
    [0.442374288, 0.141554714, 0.014211858],
    [0.491249232, 0.127914394, 0.019294232],
    [0.539365767, 0.051472401, 0.030840447]
]

enpichi2FitsAULoffset_Fa22 = [
    [0.113905344, 0.138968565, 0.119472128],
    [0.152729530, -0.029264192, 0.022574927],
    [0.198321988, -0.027342472, 0.012204149],
    [0.246523461, -0.022473325, 0.009491199],
    [0.295697442, -0.010038209, 0.008781214],
    [0.344509236, -0.014477763, 0.008665829],
    [0.393376160, -0.011451293, 0.009979557],
    [0.442374288, -0.023287795, 0.012612311],
    [0.491249232, -0.010266991, 0.017149349],
    [0.539365767, -0.021790305, 0.029119405]
]

enpichi2FitsAULsinphi_Fa22 = [
    [0.113905344, -0.307365471, 0.107178323],
    [0.152729530, -0.029496984, 0.019269801],
    [0.198321988, 0.003151551, 0.010048977],
    [0.246523461, 0.046920780, 0.007755757],
    [0.295697442, 0.071098352, 0.007218677],
    [0.344509236, 0.091948355, 0.007184757],
    [0.393376160, 0.104848180, 0.008058470],
    [0.442374288, 0.118107002, 0.010252536],
    [0.491249232, 0.061793537, 0.013639237],
    [0.539365767, 0.038433042, 0.023577991]
]

enpichi2FitsAULsin2phi_Fa22 = [
    [0.113905344, -0.760387208, 0.280519422],
    [0.152729530, -0.099007856, 0.044154736],
    [0.198321988, -0.072659147, 0.022092156],
    [0.246523461, -0.099833549, 0.016498076],
    [0.295697442, -0.134810196, 0.015230285],
    [0.344509236, -0.084765129, 0.014501825],
    [0.393376160, -0.064714516, 0.016521397],
    [0.442374288, -0.079632538, 0.020869120],
    [0.491249232, -0.074661412, 0.028755318],
    [0.539365767, -0.087239803, 0.049628509]
]

enpichi2FitsALL_Fa22 = [
    [0.113905344, 0.594624126, 0.197442658],
    [0.152729530, 0.170932120, 0.040615602],
    [0.198321988, 0.357131203, 0.027273435],
    [0.246523461, 0.418501277, 0.024797137],
    [0.295697442, 0.543063485, 0.026845723],
    [0.344509236, 0.649974056, 0.029787877],
    [0.393376160, 0.658934623, 0.032087217],
    [0.442374288, 0.727209167, 0.037947536],
    [0.491249232, 0.843410755, 0.048898255],
    [0.539365767, 0.701668739, 0.068231907]
]

enpichi2FitsALLcosphi_Fa22 = [
    [0.113905344, 2.465693219, 0.340905938],
    [0.152729530, 0.245275240, 0.066925718],
    [0.198321988, 0.268932201, 0.041933554],
    [0.246523461, 0.144995775, 0.037407799],
    [0.295697442, 0.116390203, 0.039502357],
    [0.344509236, 0.002661954, 0.043573722],
    [0.393376160, -0.040330615, 0.048033863],
    [0.442374288, -0.081585401, 0.055752838],
    [0.491249232, -0.089079540, 0.072980403],
    [0.539365767, 0.003321235, 0.101964724]
]

# RGC Sp23 data (latest fit)
enpichi2FitsALUoffset_Sp23 = [
    [0.108719696, -0.024631772, 0.013597918],
    [0.151276237, -0.008141918, 0.009183845],
    [0.198001764, 0.001754568, 0.005408275],
    [0.245456398, 0.007515013, 0.004344560],
    [0.294556806, -0.007908606, 0.004432656],
    [0.344078030, -0.003446401, 0.004791130],
    [0.393097899, 0.013336555, 0.005682381],
    [0.442440193, -0.001226365, 0.007375181],
    [0.491348572, -0.015587450, 0.010598170],
    [0.539609734, 0.002479266, 0.018419385]
]

enpichi2FitsALUsinphi_Sp23 = [
    [0.108719696, 2.183765525, 0.032065459],
    [0.151276237, 0.091076155, 0.023440237],
    [0.198001764, 0.108530693, 0.015753699],
    [0.245456398, 0.102746081, 0.013984734],
    [0.294556806, 0.140383518, 0.014322018],
    [0.344078030, 0.125626942, 0.015711544],
    [0.393097899, 0.108330407, 0.018334063],
    [0.442440193, 0.114800677, 0.022845117],
    [0.491348572, 0.099835921, 0.030471119],
    [0.539609734, 0.105807345, 0.050823998]
]

enpichi2FitsAULoffset_Sp23 = [
    [0.108719696, -0.053701744, 0.079640563],
    [0.151276237, -0.106423133, 0.026817098],
    [0.198001764, -0.033559759, 0.015578736],
    [0.245456398, -0.045808313, 0.012933798],
    [0.294556806, -0.003314873, 0.012805201],
    [0.344078030, 0.001442509, 0.013887418],
    [0.393097899, 0.019514895, 0.016651146],
    [0.442440193, 0.043299274, 0.020929924],
    [0.491348572, 0.048234170, 0.030218629],
    [0.539609734, 0.040513659, 0.051488917]
]

enpichi2FitsAULsinphi_Sp23 = [
    [0.108719696, 0.356133470, 0.084895439],
    [0.151276237, -0.068862718, 0.021654395],
    [0.198001764, -0.002620200, 0.011583065],
    [0.245456398, 0.030993754, 0.009440343],
    [0.294556806, 0.032459127, 0.009686597],
    [0.344078030, 0.053537944, 0.010656723],
    [0.393097899, 0.081146769, 0.013210087],
    [0.442440193, 0.054906223, 0.016386488],
    [0.491348572, 0.082340744, 0.023893347],
    [0.539609734, 0.006405192, 0.041683016]
]

enpichi2FitsAULsin2phi_Sp23 = [
    [0.108719696, 1.128354843, 0.185230187],
    [0.151276237, -0.067198279, 0.049009101],
    [0.198001764, -0.132055138, 0.025596539],
    [0.245456398, -0.145086096, 0.020712630],
    [0.294556806, -0.115896163, 0.020692891],
    [0.344078030, -0.109998036, 0.022059260],
    [0.393097899, -0.106358168, 0.026890976],
    [0.442440193, -0.074389542, 0.033772083],
    [0.491348572, -0.034211728, 0.050184187],
    [0.539609734, -0.067278296, 0.086342026]
]

enpichi2FitsALL_Sp23 = [
    [0.108719696, 1.037563694, 0.224337868],
    [0.151276237, 0.232153084, 0.069271855],
    [0.198001764, 0.282314016, 0.041199035],
    [0.245456398, 0.400024947, 0.037874262],
    [0.294556806, 0.483427766, 0.039252434],
    [0.344078030, 0.633175197, 0.045243568],
    [0.393097899, 0.614636678, 0.049154451],
    [0.442440193, 0.569390188, 0.057607541],
    [0.491348572, 0.705841418, 0.075773158],
    [0.539609734, 0.835777258, 0.119691021]
]

enpichi2FitsALLcosphi_Sp23 = [
    [0.108719696, 1.532108564, 0.342284354],
    [0.151276237, 0.504702592, 0.112784933],
    [0.198001764, 0.150929466, 0.065331214],
    [0.245456398, 0.188738748, 0.058328956],
    [0.294556806, 0.053344407, 0.061499368],
    [0.344078030, -0.033208604, 0.069503691],
    [0.393097899, -0.040407681, 0.072738899],
    [0.442440193, -0.080940836, 0.085208090],
    [0.491348572, -0.178354853, 0.115199973],
    [0.539609734, -0.349541366, 0.177300356]
]

# -----------------------------------------------------------------------------
# Organize data by run period
# -----------------------------------------------------------------------------
periods = {
    "Su22": {
        "ALUsinphi": enpichi2FitsALUsinphi_Su22,
        "AULsinphi": enpichi2FitsAULsinphi_Su22,
        "AULsin2phi": enpichi2FitsAULsin2phi_Su22,
        "ALL_n0": enpichi2FitsALL_Su22,
        "ALLcosphi": enpichi2FitsALLcosphi_Su22
    },
    "Fa22": {
        "ALUsinphi": enpichi2FitsALUsinphi_Fa22,
        "AULsinphi": enpichi2FitsAULsinphi_Fa22,
        "AULsin2phi": enpichi2FitsAULsin2phi_Fa22,
        "ALL_n0": enpichi2FitsALL_Fa22,
        "ALLcosphi": enpichi2FitsALLcosphi_Fa22
    },
    "Sp23": {
        "ALUsinphi": enpichi2FitsALUsinphi_Sp23,
        "AULsinphi": enpichi2FitsAULsinphi_Sp23,
        "AULsin2phi": enpichi2FitsAULsin2phi_Sp23,
        "ALL_n0": enpichi2FitsALL_Sp23,
        "ALLcosphi": enpichi2FitsALLcosphi_Sp23
    }
}

# Colors for each run period
colors = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green"
}

# -----------------------------------------------------------------------------
# Convert each list in each period to NumPy arrays
# -----------------------------------------------------------------------------
for p in periods:
    for key in periods[p]:
        data_list = periods[p][key]
        if len(data_list) > 0:
            arr = np.array(data_list)
            x = arr[:, 0]
            y = arr[:, 1]
            yerr = arr[:, 2]
            periods[p][key] = {"x": x, "y": y, "yerr": yerr}
        else:
            periods[p][key] = None
        # endif
    # end for key
# end for p

# -----------------------------------------------------------------------------
# Create output directory if it does not exist
# -----------------------------------------------------------------------------
out_dir = os.path.join("output", "enpi+")
if not os.path.isdir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
# endif

# -----------------------------------------------------------------------------
# Plotting: 1×3 figure, all three run periods on each subplot
# -----------------------------------------------------------------------------
plt.figure(figsize=(15, 5))
plt.suptitle(
    r"$ep \rightarrow en\pi^{+}$, $0.07 < |t| < 0.7$, $z > 0.55$, $y < 0.65$, $0.75 < M_{x}^{2} < 1.05\ (\mathrm{GeV}^{2})$",
    fontsize=16,
    y=0.96
)

# Increase base font size for axes labels
label_fontsize = 13

# -------------------------
# Subplot 1: ALU sinφ (all 3 periods)
# -------------------------
ax1 = plt.subplot(1, 3, 1)

for p in ["Su22", "Fa22", "Sp23"]:
    data = periods[p]["ALUsinphi"]
    if data is not None:
        ax1.errorbar(
            data["x"],
            data["y"],
            yerr=data["yerr"],
            fmt="o",
            color=colors[p],
            ecolor=colors[p],
            capsize=3,
            label=p
        )
    # endif
# end for

ax1.set_xlim(0, 0.7)
ax1.set_ylim(-0.2, 0.2)
ax1.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax1.set_ylabel(r"$F_{LU}^{\sin\phi}/F_{UU}$", fontsize=label_fontsize)
ax1.axhline(0, color="black", linestyle="--", linewidth=1.2)
ax1.grid(True, linestyle="--", alpha=0.6)

# Legend for run periods (Su22, Fa22, Sp23)
legend1 = ax1.legend(
    title="Run Period",
    frameon=True,
    edgecolor="black",
    fontsize=11,
    title_fontsize=12
)
legend1.get_frame().set_alpha(0.9)

# -------------------------
# Subplot 2: AUL sinφ (n=1, open) & sin2φ (n=2, closed) for all periods
# -------------------------
ax2 = plt.subplot(1, 3, 2, sharex=ax1)

for p in ["Su22", "Fa22", "Sp23"]:
    d1 = periods[p]["AULsinphi"]
    d2 = periods[p]["AULsin2phi"]
    if d1 is not None:
        ax2.errorbar(
            d1["x"],
            d1["y"],
            yerr=d1["yerr"],
            fmt="o",
            mfc="none",        # open circle for n=1
            mec=colors[p],
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=1"
        )
    # endif
    if d2 is not None:
        ax2.errorbar(
            d2["x"],
            d2["y"],
            yerr=d2["yerr"],
            fmt="o",
            color=colors[p],   # filled circle for n=2
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=2"
        )
    # endif
# end for

ax2.set_xlim(0, 0.7)
ax2.set_ylim(-0.2, 0.2)
ax2.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax2.set_ylabel(r"$F_{UL}^{\sin\,n\phi}/F_{UU}$", fontsize=label_fontsize)
ax2.axhline(0, color="black", linestyle="--", linewidth=1.2)
ax2.grid(True, linestyle="--", alpha=0.6)

# Legend for harmonic n (subplot 2)
legend_n2 = ax2.legend(
    handles=[
        Line2D([0], [0], marker='o', mfc='none', mec='black', linestyle='', label='n=1'),
        Line2D([0], [0], marker='o', color='black', linestyle='', label='n=2')
    ],
    title="Harmonic n",
    frameon=True,
    edgecolor="black",
    loc='upper right',
    fontsize=11,
    title_fontsize=12
)
legend_n2.get_frame().set_alpha(0.9)
ax2.add_artist(legend_n2)

# Legend for run periods (subplot 2)
legend_runs2 = ax2.legend(
    handles=[
        Line2D([0], [0], marker='o', color=colors["Su22"], linestyle='', label='Su22'),
        Line2D([0], [0], marker='o', color=colors["Fa22"], linestyle='', label='Fa22'),
        Line2D([0], [0], marker='o', color=colors["Sp23"], linestyle='', label='Sp23')
    ],
    title="Run Period",
    frameon=True,
    edgecolor="black",
    loc='lower right',
    fontsize=11,
    title_fontsize=12
)
legend_runs2.get_frame().set_alpha(0.9)

# -------------------------
# Subplot 3: ALL n=0 (open) & cosφ (n=1, closed) for all periods
# -------------------------
ax3 = plt.subplot(1, 3, 3, sharex=ax1)

for p in ["Su22", "Fa22", "Sp23"]:
    d0 = periods[p]["ALL_n0"]
    d1 = periods[p]["ALLcosphi"]
    if d0 is not None:
        ax3.errorbar(
            d0["x"],
            d0["y"],
            yerr=d0["yerr"],
            fmt="o",
            mfc="none",        # open circle for n=0
            mec=colors[p],
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=0"
        )
    # endif
    if d1 is not None:
        ax3.errorbar(
            d1["x"],
            d1["y"],
            yerr=d1["yerr"],
            fmt="o",
            color=colors[p],   # filled circle for n=1
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=1"
        )
    # endif
# end for

ax3.set_xlim(0, 0.7)
ax3.set_ylim(-1, 1)
ax3.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax3.set_ylabel(r"$F_{LL}^{\cos\,n\phi}/F_{UU}$", fontsize=label_fontsize)
ax3.axhline(0, color="black", linestyle="--", linewidth=1.2)
ax3.grid(True, linestyle="--", alpha=0.6)

# Legend for harmonic n (subplot 3)
legend_n3 = ax3.legend(
    handles=[
        Line2D([0], [0], marker='o', mfc='none', mec='black', linestyle='', label='n=0'),
        Line2D([0], [0], marker='o', color='black', linestyle='', label='n=1')
    ],
    title="Harmonic n",
    frameon=True,
    edgecolor="black",
    loc='upper right',
    fontsize=11,
    title_fontsize=12
)
legend_n3.get_frame().set_alpha(0.9)
ax3.add_artist(legend_n3)

# Legend for run periods (subplot 3)
legend_runs3 = ax3.legend(
    handles=[
        Line2D([0], [0], marker='o', color=colors["Su22"], linestyle='', label='Su22'),
        Line2D([0], [0], marker='o', color=colors["Fa22"], linestyle='', label='Fa22'),
        Line2D([0], [0], marker='o', color=colors["Sp23"], linestyle='', label='Sp23')
    ],
    title="Run Period",
    frameon=True,
    edgecolor="black",
    loc='lower right',
    fontsize=11,
    title_fontsize=12
)
legend_runs3.get_frame().set_alpha(0.9)

plt.tight_layout(rect=[0, 0, 1, 0.93])

# Save the asymmetry figure under output/enpi+/
asymmetry_filename = os.path.join(out_dir, "rgc_enpi+_AllPeriods.pdf")
plt.savefig(asymmetry_filename)
print(f"Asymmetry plot saved to '{asymmetry_filename}'")

# -----------------------------------------------------------------------------
# New functionality: Plot dilution factor vs xB for Su22 and Fa22
# -----------------------------------------------------------------------------

# x–values from the Su22 sinφ asymmetries (10 bins)
x_Su22 = np.array([row[0] for row in enpichi2FitsALUsinphi_Su22])

# New Su22 dilution factors (10 bins)
dil_Su22     = np.array([0.482784, 0.311356, 0.325547, 0.354333, 0.355804,
                         0.362807, 0.374695, 0.394802, 0.331315, 0.398176])
dil_err_Su22 = np.array([0.11728,  0.0380967, 0.0195395, 0.0139435, 0.0118704,
                         0.0113884, 0.012495,  0.015489,  0.0263908, 0.0381037])

# x–values from the Fa22 sinφ asymmetries (10 bins)
x_Fa22 = np.array([row[0] for row in enpichi2FitsALUsinphi_Fa22])

# New Fa22 dilution factors (10 bins)
dil_Fa22     = np.array([0.387857, 0.360117, 0.348367, 0.3514,   0.343512,
                         0.360489, 0.353206, 0.366436, 0.343407, 0.358797])
dil_err_Fa22 = np.array([0.056521, 0.0125693, 0.00699825, 0.00527699, 0.00469843,
                         0.0043642,  0.00512107, 0.00646144, 0.0100004,  0.0164314])

plt.figure(figsize=(6, 5))
plt.errorbar(
    x_Su22, dil_Su22, yerr=dil_err_Su22,
    fmt="o", color=colors["Su22"], ecolor=colors["Su22"],
    capsize=3, label="Su22"
)
plt.errorbar(
    x_Fa22, dil_Fa22, yerr=dil_err_Fa22,
    fmt="o", color=colors["Fa22"], ecolor=colors["Fa22"],
    capsize=3, label="Fa22"
)

plt.xlabel(r"$x_{B}$", fontsize=label_fontsize)
plt.ylabel(r"$D_{f}$", fontsize=label_fontsize)
plt.ylim(0, 1)
plt.xlim(0, 0.7)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(frameon=True, edgecolor="black",
           fontsize=11, title="Run Period", title_fontsize=12)

plt.tight_layout()
dilution_filename = os.path.join(out_dir, "dilution_vs_xB.pdf")
plt.savefig(dilution_filename)
print(f"Dilution factor plot saved to '{dilution_filename}'")

# -----------------------------------------------------------------------------
# New functionality: Plot mean Q², W, z, and x_F vs x_B for Su22, Fa22, and Sp23
# -----------------------------------------------------------------------------

# Mean kinematic values for Su22
xB_Su22 = np.array([0.094, 0.169, 0.255, 0.348, 0.441, 0.535])
Q2_Su22 = np.array([1.325, 1.862, 2.240, 2.592, 3.374, 4.589])
W_Su22  = np.array([3.689, 3.179, 2.722, 2.381, 2.256, 2.197])
z_Su22  = np.array([0.972, 0.971, 0.955, 0.932, 0.919, 0.912])
xF_Su22 = np.array([0.883, 0.862, 0.812, 0.756, 0.745, 0.757])

# Mean kinematic values for Fa22
xB_Fa22 = np.array([0.095, 0.168, 0.255, 0.348, 0.441, 0.535])
Q2_Fa22 = np.array([1.341, 1.874, 2.254, 2.612, 3.388, 4.607])
W_Fa22  = np.array([3.698, 3.193, 2.728, 2.389, 2.261, 2.200])
z_Fa22  = np.array([0.973, 0.969, 0.953, 0.930, 0.918, 0.912])
xF_Fa22 = np.array([0.886, 0.861, 0.810, 0.755, 0.744, 0.756])

# Mean kinematic values for Sp23
xB_Sp23 = np.array([0.091, 0.167, 0.251, 0.346, 0.441, 0.535])
Q2_Sp23 = np.array([1.132, 1.502, 1.864, 2.513, 3.379, 4.607])
W_Sp23  = np.array([3.487, 2.896, 2.517, 2.356, 2.258, 2.199])
z_Sp23  = np.array([0.961, 0.960, 0.943, 0.929, 0.919, 0.912])
xF_Sp23 = np.array([0.854, 0.824, 0.771, 0.750, 0.745, 0.757])

# Create a 2×2 figure
fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# Common plot settings
marker_Su22 = 'o'
marker_Fa22 = '^'
marker_Sp23 = 's'
color_Su22 = colors["Su22"]
color_Fa22 = colors["Fa22"]
color_Sp23 = colors["Sp23"]
label_fontsize = 12
tick_fontsize = 10

# Top-left: Q² vs x_B
ax = axes[0, 0]
ax.plot(xB_Su22, Q2_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='', label='Su22')
ax.plot(xB_Fa22, Q2_Fa22, marker_Fa22, color=color_Fa22, markersize=6, linestyle='', label='Fa22')
ax.plot(xB_Sp23, Q2_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle Q^{2} \rangle\ \mathrm{(GeV^{2})}$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Top-right: W vs x_B
ax = axes[0, 1]
ax.plot(xB_Su22, W_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='', label='Su22')
ax.plot(xB_Fa22, W_Fa22, marker_Fa22, color=color_Fa22, markersize=6, linestyle='', label='Fa22')
ax.plot(xB_Sp23, W_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle W \rangle\ \mathrm{(GeV)}$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Bottom-left: z vs x_B
ax = axes[1, 0]
ax.plot(xB_Su22, z_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='', label='Su22')
ax.plot(xB_Fa22, z_Fa22, marker_Fa22, color=color_Fa22, markersize=6, linestyle='', label='Fa22')
ax.plot(xB_Sp23, z_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle z \rangle$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.set_ylim(0.90, 1.0)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Bottom-right: x_F vs x_B
ax = axes[1, 1]
ax.plot(xB_Su22, xF_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='', label='Su22')
ax.plot(xB_Fa22, xF_Fa22, marker_Fa22, color=color_Fa22, markersize=6, linestyle='', label='Fa22')
ax.plot(xB_Sp23, xF_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle x_{F} \rangle$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Add a single legend for the entire figure
handles = [
    Line2D([0], [0], marker=marker_Su22, color=color_Su22, linestyle='', label='Su22', markersize=6),
    Line2D([0], [0], marker=marker_Fa22, color=color_Fa22, linestyle='', label='Fa22', markersize=6),
    Line2D([0], [0], marker=marker_Sp23, color=color_Sp23, linestyle='', label='Sp23', markersize=6)
]
fig.legend(
    handles=handles,
    loc='upper center',
    ncol=3,
    frameon=True,
    edgecolor="black",
    fontsize=11,
    title="Data Set",
    title_fontsize=12
)

plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save the kinematic comparison figure
kinematic_filename = os.path.join(out_dir, "kinematic_comparison_Su22_Fa22_Sp23.pdf")
plt.savefig(kinematic_filename)
print(f"Kinematic comparison plot saved as '{kinematic_filename}'")


#########
def compare_fa22_versions():
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # New, high-contrast colors
    cmp_colors = {"on": "tab:purple", "off": "tab:cyan"}

    # The two Fa22 configurations
    fa22_on = {
        "ALUsinphi": enpichi2FitsALUsinphi_Fa22,
        "AULsinphi": enpichi2FitsAULsinphi_Fa22,
        "AULsin2phi": enpichi2FitsAULsin2phi_Fa22,
        "ALL_n0":     enpichi2FitsALL_Fa22,
        "ALLcosphi":  enpichi2FitsALLcosphi_Fa22
    }
    fa22_off = {
        "ALUsinphi": [
            [0.113810624, 0.247981539, 0.100196499],
            [0.152990304, 0.082747546, 0.019644028],
            [0.198503925, 0.105733843, 0.011966558],
            [0.246596319, 0.126546450, 0.010121316],
            [0.295736697, 0.139206646, 0.009854151],
            [0.344508054, 0.140353985, 0.010319784],
            [0.393387120, 0.141851722, 0.011648971],
            [0.442375088, 0.144403863, 0.014456324],
            [0.491252516, 0.126450678, 0.019621998],
            [0.539338687, 0.058299070, 0.031464698]
        ],
        "AULsinphi": [
            [0.113810624, -0.173630451, 0.135281218],
            [0.152990304, -0.041016757, 0.022578625],
            [0.198503925, 0.014087512, 0.011434166],
            [0.246596319, 0.058831967, 0.009079952],
            [0.295736697, 0.087451423, 0.008763746],
            [0.344508054, 0.105809831, 0.008747582],
            [0.393387120, 0.123913717, 0.009818065],
            [0.442375088, 0.136204431, 0.012202708],
            [0.491252516, 0.071706367, 0.015237049],
            [0.539338687, 0.044687045, 0.026092727]
        ],
        "AULsin2phi": [
            [0.113810624, -0.176545812, 0.377457035],
            [0.152990304, -0.173650111, 0.052737531],
            [0.198503925, -0.076016447, 0.025356810],
            [0.246596319, -0.103297299, 0.019149828],
            [0.295736697, -0.143183202, 0.018477227],
            [0.344508054, -0.091521832, 0.017431681],
            [0.393387120, -0.066784555, 0.020019969],
            [0.442375088, -0.083773315, 0.024427119],
            [0.491252516, -0.078812199, 0.032337456],
            [0.539338687, -0.085419563, 0.054973279]
        ],
        "ALL_n0": [
            [0.113810624, -0.028578533, 0.296795894],
            [0.152990304, 0.177588470, 0.048784416],
            [0.198503925, 0.400180058, 0.033003491],
            [0.246596319, 0.453029971, 0.030200450],
            [0.295736697, 0.591420857, 0.033664672],
            [0.344508054, 0.708929957, 0.037813759],
            [0.393387120, 0.716648913, 0.040078956],
            [0.442375088, 0.789795414, 0.046697140],
            [0.491252516, 0.897801126, 0.058778724],
            [0.539338687, 0.777727415, 0.078360736]
        ],
        "ALLcosphi": [
            [0.113810624, 1.192146786, 0.553554921],
            [0.152990304, 0.279249263, 0.081273811],
            [0.198503925, 0.300263142, 0.050320471],
            [0.246596319, 0.171827757, 0.045498535],
            [0.295736697, 0.141019633, 0.049337298],
            [0.344508054, -0.012470967, 0.055330505],
            [0.393387120, -0.044620314, 0.060446715],
            [0.442375088, -0.070996241, 0.068981120],
            [0.491252516, -0.132769930, 0.087742825],
            [0.539338687, 0.030052764, 0.117432841]
        ]
    }

    # helper: list-of-lists → dict with x,y,yerr
    def to_dict(data):
        a = np.array(data)
        return {"x": a[:,0], "y": a[:,1], "yerr": a[:,2]}

    on  = {k: to_dict(v) for k,v in fa22_on.items()}
    off = {k: to_dict(v) for k,v in fa22_off.items()}

    fig, axes = plt.subplots(1, 3, figsize=(15,5))
    plt.suptitle("Fa22: Fiducial Cuts ON vs OFF", fontsize=16, y=0.97)
    fs = 13

    # --- 1) ALU sinφ ---
    ax = axes[0]
    ax.errorbar(on["ALUsinphi"]["x"], on["ALUsinphi"]["y"],  yerr=on["ALUsinphi"]["yerr"],
                fmt="o", color=cmp_colors["on"],  capsize=3, label="cuts ON")
    ax.errorbar(off["ALUsinphi"]["x"], off["ALUsinphi"]["y"], yerr=off["ALUsinphi"]["yerr"],
                fmt="o", color=cmp_colors["off"], capsize=3, label="cuts OFF")
    ax.set(xlim=(0,0.7), ylim=(-0.2,0.2),
           xlabel=r"$x_B$", ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    ax.axhline(0, linestyle="--", color="black")
    ax.grid(True,linestyle="--",alpha=0.5)
    ax.legend(title="Configuration", frameon=True)

    # --- 2) AUL sinφ (open) & sin2φ (filled) ---
    ax = axes[1]
    # ON
    ax.errorbar(on["AULsinphi"]["x"], on["AULsinphi"]["y"], yerr=on["AULsinphi"]["yerr"],
                fmt="o", mfc="none", mec=cmp_colors["on"], capsize=3)
    ax.errorbar(on["AULsin2phi"]["x"],on["AULsin2phi"]["y"],yerr=on["AULsin2phi"]["yerr"],
                fmt="o", color=cmp_colors["on"], capsize=3)
    # OFF
    ax.errorbar(off["AULsinphi"]["x"], off["AULsinphi"]["y"], yerr=off["AULsinphi"]["yerr"],
                fmt="s", mfc="none", mec=cmp_colors["off"], capsize=3)
    ax.errorbar(off["AULsin2phi"]["x"],off["AULsin2phi"]["y"],yerr=off["AULsin2phi"]["yerr"],
                fmt="s", color=cmp_colors["off"], capsize=3)

    ax.set(xlim=(0,0.7), ylim=(-0.2,0.2),
           xlabel=r"$x_B$", ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    ax.axhline(0, linestyle="--", color="black")
    ax.grid(True,linestyle="--",alpha=0.5)
    # harmonic legend
    harm = ax.legend(
        handles=[
            Line2D([0],[0], marker='o', mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker='o', color='black', linestyle='', label='n=2')
        ],
        title="Harmonic", loc='upper right', frameon=True
    )
    ax.add_artist(harm)
    # config legend
    ax.legend(
        handles=[
            Line2D([0],[0], marker='o', color=cmp_colors["on"], linestyle='', label='cuts ON'),
            Line2D([0],[0], marker='s', color=cmp_colors["off"], linestyle='', label='cuts OFF')
        ],
        title="Config", loc='lower right', frameon=True
    )

    # --- 3) ALL n=0 (open) & cosφ (filled) ---
    ax = axes[2]
    # ON
    ax.errorbar(on["ALL_n0"]["x"], on["ALL_n0"]["y"], yerr=on["ALL_n0"]["yerr"],
                fmt="o", mfc="none", mec=cmp_colors["on"], capsize=3)
    ax.errorbar(on["ALLcosphi"]["x"], on["ALLcosphi"]["y"], yerr=on["ALLcosphi"]["yerr"],
                fmt="o", color=cmp_colors["on"], capsize=3)
    # OFF
    ax.errorbar(off["ALL_n0"]["x"], off["ALL_n0"]["y"], yerr=off["ALL_n0"]["yerr"],
                fmt="s", mfc="none", mec=cmp_colors["off"], capsize=3)
    ax.errorbar(off["ALLcosphi"]["x"], off["ALLcosphi"]["y"], yerr=off["ALLcosphi"]["yerr"],
                fmt="s", color=cmp_colors["off"], capsize=3)

    ax.set(xlim=(0,0.7), ylim=(-1,1),
           xlabel=r"$x_B$", ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    ax.axhline(0, linestyle="--", color="black")
    ax.grid(True,linestyle="--",alpha=0.5)
    # harmonic legend
    harm3 = ax.legend(
        handles=[
            Line2D([0],[0], marker='o', mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker='o', color='black', linestyle='', label='n=1')
        ],
        title="Harmonic", loc='upper right', frameon=True
    )
    ax.add_artist(harm3)
    # config legend
    ax.legend(
        handles=[
            Line2D([0],[0], marker='o', color=cmp_colors["on"], linestyle='', label='cuts ON'),
            Line2D([0],[0], marker='s', color=cmp_colors["off"], linestyle='', label='cuts OFF')
        ],
        title="Config", loc='lower right', frameon=True
    )

    plt.tight_layout(rect=[0,0,1,0.94])

    # <-- Save before showing -->
    cmp_file = os.path.join(out_dir, "fa22_cuts_on_vs_off.pdf")
    plt.savefig(cmp_file)
    print(f"Saved comparison plot to '{cmp_file}'")

    plt.show()

# Call it at the end
compare_fa22_versions()