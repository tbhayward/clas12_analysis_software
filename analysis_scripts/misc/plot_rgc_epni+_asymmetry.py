#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Hard-coded data for all three run periods
#
# UPDATED Su22, Fa22, Sp23 asymmetry results (10 bins each)
# Each row: [xB, value, stat_err]
# -----------------------------------------------------------------------------

# RGC Su22 data (re-analysis)
enpichi2FitsALUoffset_Su22 = [
    [0.110552601, -0.006117919, 0.000215170],
    [0.150741944,  0.000078444, 0.000110128],
    [0.197673274, -0.000003040, 0.000069016],
    [0.246387600,  0.000135649, 0.000054752],
    [0.295782734,  0.000104310, 0.000048573],
    [0.344980771,  0.000030123, 0.000046113],
    [0.393720730,  0.000079925, 0.000050308],
    [0.443065481,  0.000021758, 0.000061083],
    [0.492095071,  0.000029211, 0.000081769],
    [0.541127567, -0.000037368, 0.000119702]
]

enpichi2FitsALUsinphi_Su22 = [
    [0.110552601, 0.012057485, 0.000438439],
    [0.150741944, 0.000906295, 0.000223347],
    [0.197673274, 0.000830652, 0.000158988],
    [0.246387600, 0.001082510, 0.000140561],
    [0.295782734, 0.001503434, 0.000138592],
    [0.344980771, 0.001516246, 0.000141378],
    [0.393720730, 0.001356702, 0.000156362],
    [0.443065481, 0.001171826, 0.000183889],
    [0.492095071, 0.001149614, 0.000232640],
    [0.541127567, 0.000979680, 0.000316091]
]

enpichi2FitsAULoffset_Su22 = [
    [0.110552601, 0.044741243, 0.099349046],
    [0.150741944, 0.062425143, 0.029925156],
    [0.197673274, 0.061671858, 0.018718475],
    [0.246387600, 0.034964256, 0.014584925],
    [0.295782734, 0.036894781, 0.013475149],
    [0.344980771, 0.030965579, 0.013041849],
    [0.393720730, 0.023976428, 0.013893082],
    [0.443065481, 0.011776591, 0.016547362],
    [0.492095071, 0.034442169, 0.022061922],
    [0.541127567, 0.029570685, 0.031112595]
]

enpichi2FitsAULsinphi_Su22 = [
    [0.110552601, -0.018648666, 0.109429566],
    [0.150741944, -0.084714937, 0.027493695],
    [0.197673274, -0.054382422, 0.015892852],
    [0.246387600, -0.015151252, 0.011794335],
    [0.295782734,  0.021777436, 0.010463599],
    [0.344980771,  0.039896688, 0.010366710],
    [0.393720730,  0.059814159, 0.011178480],
    [0.443065481,  0.066655670, 0.013666322],
    [0.492095071,  0.082821514, 0.017892727],
    [0.541127567,  0.025017433, 0.025137572]
]

enpichi2FitsAULsin2phi_Su22 = [
    [0.110552601,  0.453959312, 0.284292831],
    [0.150741944, -0.046225923, 0.061845582],
    [0.197673274, -0.053177463, 0.034452445],
    [0.246387600, -0.094313020, 0.025628178],
    [0.295782734, -0.129174732, 0.023375259],
    [0.344980771, -0.128732206, 0.022069033],
    [0.393720730, -0.072701977, 0.023847765],
    [0.443065481, -0.093730861, 0.028216138],
    [0.492095071, -0.093554426, 0.037162619],
    [0.541127567, -0.106219441, 0.052629739]
]

enpichi2FitsALL_Su22 = [
    [0.110552601, 0.002509421, 0.001915880],
    [0.150741944, 0.002403557, 0.000524823],
    [0.197673274, 0.003365131, 0.000399786],
    [0.246387600, 0.004338824, 0.000402580],
    [0.295782734, 0.005188224, 0.000439416],
    [0.344980771, 0.005782536, 0.000475841],
    [0.393720730, 0.006725732, 0.000540518],
    [0.443065481, 0.007886588, 0.000639369],
    [0.492095071, 0.008983387, 0.000751583],
    [0.541127567, 0.009581715, 0.000893905]
]

enpichi2FitsALLcosphi_Su22 = [
    [0.110552601,  0.008106012, 0.003335305],
    [0.150741944,  0.003264510, 0.000900516],
    [0.197673274,  0.002402400, 0.000619791],
    [0.246387600,  0.001139108, 0.000614646],
    [0.295782734,  0.000541307, 0.000648233],
    [0.344980771, -0.000219382, 0.000751314],
    [0.393720730, -0.000143300, 0.000859220],
    [0.443065481,  0.000734571, 0.000969747],
    [0.492095071,  0.000262492, 0.001151187],
    [0.541127567, -0.000364386, 0.001391548]
]

# RGC Fa22 data (re-analysis)
enpichi2FitsALUoffset_Fa22 = [
    [0.110589889, -0.010374073, 0.000091105],
    [0.150451165,  0.000034352, 0.000088034],
    [0.197681250, -0.000022418, 0.000055866],
    [0.246401288, -0.000045283, 0.000044262],
    [0.295848133, -0.000042643, 0.000039364],
    [0.344998460,  0.000060688, 0.000037529],
    [0.393879970,  0.000043022, 0.000040176],
    [0.442866756,  0.000003981, 0.000048849],
    [0.492158528, -0.000005131, 0.000064980],
    [0.541133012,  0.000004333, 0.000095002]
]

enpichi2FitsALUsinphi_Fa22 = [
    [0.110589889, 0.008719613, 0.000505830],
    [0.150451165, 0.000687857, 0.000176757],
    [0.197681250, 0.000865779, 0.000127211],
    [0.246401288, 0.001147578, 0.000113361],
    [0.295848133, 0.001329243, 0.000110770],
    [0.344998460, 0.001337965, 0.000113707],
    [0.393879970, 0.001301079, 0.000123840],
    [0.442866756, 0.001417225, 0.000146265],
    [0.492158528, 0.001373436, 0.000183556],
    [0.541133012, 0.001193884, 0.000251317]
]

enpichi2FitsAULoffset_Fa22 = [
    [0.110589889, -0.196012647, 0.065470550],
    [0.150451165,  0.023255806, 0.020641754],
    [0.197681250,  0.017193250, 0.013063807],
    [0.246401288,  0.001342416, 0.010356460],
    [0.295848133,  0.013080991, 0.009334176],
    [0.344998460,  0.013413240, 0.008959808],
    [0.393879970,  0.013233664, 0.009589846],
    [0.442866756,  0.001815413, 0.011758375],
    [0.492158528,  0.021317912, 0.015043222],
    [0.541133012,  0.026606701, 0.022012580]
]

enpichi2FitsAULsinphi_Fa22 = [
    [0.110589889,  0.111919893, 0.076855059],
    [0.150451165, -0.064168308, 0.018813050],
    [0.197681250, -0.038350611, 0.011058122],
    [0.246401288,  0.008692876, 0.008376255],
    [0.295848133,  0.039143313, 0.007378066],
    [0.344998460,  0.064626155, 0.007045107],
    [0.393879970,  0.079402161, 0.007535947],
    [0.442866756,  0.105584414, 0.009342104],
    [0.492158528,  0.066118686, 0.011918436],
    [0.541133012,  0.054773925, 0.017931663]
]

enpichi2FitsAULsin2phi_Fa22 = [
    [0.110589889,  0.430169746, 0.202708531],
    [0.150451165, -0.057970318, 0.045799193],
    [0.197681250, -0.079633321, 0.024870967],
    [0.246401288, -0.120977720, 0.018540109],
    [0.295848133, -0.151200234, 0.016022676],
    [0.344998460, -0.120056741, 0.015030860],
    [0.393879970, -0.087572523, 0.015867571],
    [0.442866756, -0.125888713, 0.019780859],
    [0.492158528, -0.087670793, 0.025388480],
    [0.541133012, -0.102357247, 0.038036054]
]

enpichi2FitsALL_Fa22 = [
    [0.110589889, 0.006678034, 0.001130911],
    [0.150451165, 0.002411435, 0.000360353],
    [0.197681250, 0.003450798, 0.000257569],
    [0.246401288, 0.004039149, 0.000237394],
    [0.295848133, 0.005683464, 0.000255488],
    [0.344998460, 0.006047932, 0.000265844],
    [0.393879970, 0.006617697, 0.000289208],
    [0.442866756, 0.007776840, 0.000339897],
    [0.492158528, 0.008584083, 0.000408247],
    [0.541133012, 0.008271637, 0.000512672]
]

enpichi2FitsALLcosphi_Fa22 = [
    [0.110589889,  0.015885775, 0.001999932],
    [0.150451165,  0.002029953, 0.000622120],
    [0.197681250,  0.001934395, 0.000411680],
    [0.246401288,  0.001324995, 0.000370514],
    [0.295848133,  0.001141469, 0.000392560],
    [0.344998460, -0.000344331, 0.000407963],
    [0.393879970, -0.000199290, 0.000442303],
    [0.442866756, -0.000648806, 0.000516368],
    [0.492158528, -0.000294809, 0.000621550],
    [0.541133012, -0.000230185, 0.000782815]
]


enpichi2FitsALUoffset_Sp23 = [
    [0.111631653,  0.010060553, 0.000000635],
    [0.150676114,  0.000044156, 0.000209636],
    [0.197592277, -0.000078224, 0.000133233],
    [0.246413670,  0.000047716, 0.000103075],
    [0.295614241, -0.000044049, 0.000091232],
    [0.345048781,  0.000000882, 0.000087478],
    [0.393829589,  0.000066662, 0.000094382],
    [0.442940060,  0.000081650, 0.000113707],
    [0.492017751,  0.000037476, 0.000149818],
    [0.541045948, -0.000165876, 0.000222839]
]

enpichi2FitsALUsinphi_Sp23 = [
    [0.111631653, 0.036148355, 0.000001989],
    [0.150676114, 0.000723932, 0.000428144],
    [0.197592277, 0.000913464, 0.000304538],
    [0.246413670, 0.000707874, 0.000265147],
    [0.295614241, 0.001146149, 0.000258142],
    [0.345048781, 0.001393441, 0.000265420],
    [0.393829589, 0.001553303, 0.000291427],
    [0.442940060, 0.001095648, 0.000342149],
    [0.492017751, 0.000953804, 0.000426679],
    [0.541045948, 0.001636735, 0.000591666]
]

enpichi2FitsAULoffset_Sp23 = [
    [0.111631653,  0.136199442, 0.128787458],
    [0.150676114,  0.029279325, 0.045250670],
    [0.197592277,  0.014672032, 0.028363257],
    [0.246413670,  0.001745416, 0.021847173],
    [0.295614241,  0.008556640, 0.019498906],
    [0.345048781, -0.038068124, 0.018653605],
    [0.393829589,  0.018703342, 0.020118715],
    [0.442940060, -0.003989142, 0.024079028],
    [0.492017751,  0.026282286, 0.031724364],
    [0.541045948, -0.000292453, 0.047225285]
]

enpichi2FitsAULsinphi_Sp23 = [
    [0.111631653,  0.866766985, 0.137742486],
    [0.150676114, -0.100626654, 0.041489154],
    [0.197592277, -0.057611134, 0.023925203],
    [0.246413670, -0.023432921, 0.017705931],
    [0.295614241, -0.021128331, 0.015391484],
    [0.345048781,  0.039251506, 0.014530356],
    [0.393829589,  0.059627409, 0.015741805],
    [0.442940060,  0.065438015, 0.018919386],
    [0.492017751,  0.073487105, 0.025298837],
    [0.541045948,  0.049329476, 0.038234710]
]

enpichi2FitsAULsin2phi_Sp23 = [
    [0.111631653,  2.857194524, 0.331339458],
    [0.150676114,  0.016013938, 0.098012304],
    [0.197592277, -0.060002078, 0.052776132],
    [0.246413670, -0.113583932, 0.037899010],
    [0.295614241, -0.148661027, 0.032747644],
    [0.345048781, -0.141291824, 0.030509383],
    [0.393829589, -0.090959586, 0.032793067],
    [0.442940060, -0.082551691, 0.039070699],
    [0.492017751, -0.024470892, 0.052579604],
    [0.541045948, -0.051417857, 0.079799733]
]

enpichi2FitsALL_Sp23 = [
    [0.111631653, 0.007277864, 0.001659476],
    [0.150676114, 0.001124563, 0.000787341],
    [0.197592277, 0.002791202, 0.000541450],
    [0.246413670, 0.003855443, 0.000484003],
    [0.295614241, 0.003953090, 0.000474134],
    [0.345048781, 0.005408119, 0.000505057],
    [0.393829589, 0.006517009, 0.000558322],
    [0.442940060, 0.005696669, 0.000636848],
    [0.492017751, 0.007172942, 0.000782386],
    [0.541045948, 0.009574252, 0.001085681]
]

enpichi2FitsALLcosphi_Sp23 = [
    [0.111631653,  0.012818093, 0.002583884],
    [0.150676114,  0.003437582, 0.001339065],
    [0.197592277,  0.001734363, 0.000872856],
    [0.246413670,  0.002024529, 0.000754444],
    [0.295614241, -0.000230681, 0.000727805],
    [0.345048781, -0.000782543, 0.000765856],
    [0.393829589, -0.000238505, 0.000839981],
    [0.442940060, -0.001160708, 0.000951296],
    [0.492017751, -0.000347749, 0.001166245],
    [0.541045948, -0.000143780, 0.001645159]
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
# Plotting: 1x3 figure, all three run periods on each subplot
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
# Subplot 1: ALU sinphi (all 3 periods)
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
ax1.set_ylim(-0.02, 0.02)
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
# Subplot 2: AUL sinphi (n=1, open) and sin2phi (n=2, filled)
# -------------------------
ax2 = plt.subplot(1, 3, 2, sharex=ax1)

for p in ["Su22", "Fa22", "Sp23"]:
    d1 = periods[p]["AULsinphi"]
    d2 = periods[p]["AULsin2phi"]
    if d1 is not None:
        ax2.errorbar(
            d1["x"], d1["y"], yerr=d1["yerr"],
            fmt="o", mfc="none", mec=colors[p], ecolor=colors[p],
            capsize=3, label=f"{p}, n=1"
        )
    # endif
    if d2 is not None:
        ax2.errorbar(
            d2["x"], d2["y"], yerr=d2["yerr"],
            fmt="o", color=colors[p], ecolor=colors[p],
            capsize=3, label=f"{p}, n=2"
        )
    # endif
# end for

ax2.set_xlim(0, 0.7)
ax2.set_ylim(-0.2, 0.6)
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
# Subplot 3: ALL n=0 (open) and cosphi (n=1, filled)
# -------------------------
ax3 = plt.subplot(1, 3, 3, sharex=ax1)

for p in ["Su22", "Fa22", "Sp23"]:
    d0 = periods[p]["ALL_n0"]
    d1 = periods[p]["ALLcosphi"]
    if d0 is not None:
        ax3.errorbar(
            d0["x"], d0["y"], yerr=d0["yerr"],
            fmt="o", mfc="none", mec=colors[p], ecolor=colors[p],
            capsize=3, label=f"{p}, n=0"
        )
    # endif
    if d1 is not None:
        ax3.errorbar(
            d1["x"], d1["y"], yerr=d1["yerr"],
            fmt="o", color=colors[p], ecolor=colors[p],
            capsize=3, label=f"{p}, n=1"
        )
    # endif
# end for

ax3.set_xlim(0, 0.7)
ax3.set_ylim(-0.02, 0.02)
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
# Dilution factor plot (same aesthetics as your DF script), but integrated here
# Uses x from each period's ALUsinphi x grid (10 bins), and your new Df values
# -----------------------------------------------------------------------------

# x arrays from ALU sinphi fits (per period)
x_Su22 = periods["Su22"]["ALUsinphi"]["x"]
x_Fa22 = periods["Fa22"]["ALUsinphi"]["x"]
x_Sp23 = periods["Sp23"]["ALUsinphi"]["x"]

# Su22 dilution factors (10 bins)
dil_Su22 = np.array([0.403247, 0.425510, 0.410022, 0.395036, 0.371607,
                     0.421003, 0.419785, 0.419586, 0.381923, 0.415297])
dil_err_Su22 = np.array([0.090710, 0.0257454, 0.0180958, 0.0150382, 0.0136724,
                         0.0108513, 0.0116233, 0.0142450, 0.0227520, 0.0320867])

# Fa22 dilution factors (10 bins)
dil_Fa22 = np.array([0.429508, 0.441729, 0.422862, 0.410205, 0.410415,
                     0.419988, 0.425389, 0.432687, 0.443039, 0.451170])
dil_err_Fa22 = np.array([0.0382782, 0.0101401, 0.0069170, 0.00563764, 0.00491716,
                         0.0043941, 0.0047178, 0.0056790, 0.00759078, 0.0117991])

# Sp23 dilution factors (10 bins)
dil_Sp23 = np.array([0.470499, 0.443204, 0.413036, 0.419151, 0.412881,
                     0.423348, 0.413772, 0.443118, 0.437312, 0.475676])
dil_err_Sp23 = np.array([0.0541136, 0.0197316, 0.0129473, 0.0105135, 0.00924795,
                         0.00866588, 0.00966092, 0.0110117, 0.0152021, 0.0207904])

plt.figure(figsize=(6, 5))
plt.errorbar(
    x_Su22, dil_Su22, yerr=dil_err_Su22,
    fmt="o", color=colors["Su22"], ecolor=colors["Su22"],
    capsize=3, label="Su22"
)
plt.errorbar(
    x_Fa22, dil_Fa22, yerr=dil_err_Fa22,
    fmt="s", color=colors["Fa22"], ecolor=colors["Fa22"],
    capsize=3, label="Fa22"
)
plt.errorbar(
    x_Sp23, dil_Sp23, yerr=dil_err_Sp23,
    fmt="^", color=colors["Sp23"], ecolor=colors["Sp23"],
    capsize=3, label="Sp23"
)

plt.xlabel(r"$x_{B}$", fontsize=label_fontsize)
plt.ylabel(r"$D_{f}$", fontsize=label_fontsize)
plt.xlim(0, 0.8)
# Expanded y-range to contain new values while keeping the same look
plt.ylim(0.30, 0.50)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(frameon=True, edgecolor="black",
           fontsize=11, title="Run Period", title_fontsize=12)

plt.tight_layout()
dilution_filename = os.path.join(out_dir, "dilution_vs_xB.pdf")
plt.savefig(dilution_filename)
print(f"Dilution factor plot saved to '{dilution_filename}'")

# -----------------------------------------------------------------------------
# Mean kinematics plots (unchanged content)
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

# Create a 2x2 figure
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

# Top-left: Q2 vs xB
ax = axes[0, 0]
ax.plot(xB_Su22, Q2_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='', label='Su22')
ax.plot(xB_Fa22, Q2_Fa22, marker_Fa22, color=color_Fa22, markersize=6, linestyle='', label='Fa22')
ax.plot(xB_Sp23, Q2_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle Q^{2} \rangle\ \mathrm{(GeV^{2})}$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Top-right: W vs xB
ax = axes[0, 1]
ax.plot(xB_Su22, W_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='', label='Su22')
ax.plot(xB_Fa22, W_Fa22, marker_Fa22, color=color_Fa22, markersize=6, linestyle='', label='Fa22')
ax.plot(xB_Sp23, W_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle W \rangle\ \mathrm{(GeV)}$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Bottom-left: z vs xB
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

# Bottom-right: xF vs xB
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
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # reuse your original run-period colors
    cmp_colors = {"on": colors["Fa22"], "off": colors["Su22"]}

    # two Fa22 configurations
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
            [0.198503925,  0.014087512, 0.011434166],
            [0.246596319,  0.058831967, 0.009079952],
            [0.295736697,  0.087451423, 0.008763746],
            [0.344508054,  0.105809831, 0.008747582],
            [0.393387120,  0.123913717, 0.009818065],
            [0.442375088,  0.136204431, 0.012202708],
            [0.491252516,  0.071706367, 0.015237049],
            [0.539338687,  0.044687045, 0.026092727]
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
            [0.152990304,  0.177588470, 0.048784416],
            [0.198503925,  0.400180058, 0.033003491],
            [0.246596319,  0.453029971, 0.030200450],
            [0.295736697,  0.591420857, 0.033664672],
            [0.344508054,  0.708929957, 0.037813759],
            [0.393387120,  0.716648913, 0.040078956],
            [0.442375088,  0.789795414, 0.046697140],
            [0.491252516,  0.897801126, 0.058778724],
            [0.539338687,  0.777727415, 0.078360736]
        ],
        "ALLcosphi": [
            [0.113810624, 1.192146786, 0.553554921],
            [0.152990304, 0.279249263, 0.081273811],
            [0.198503925, 0.300263142, 0.050320471],
            [0.246596319, 0.171827757, 0.045498535],
            [0.295736697, 0.141019633, 0.049337298],
            [0.344508054,-0.012470967, 0.055330505],
            [0.393387120,-0.044620314, 0.060446715],
            [0.442375088,-0.070996241, 0.068981120],
            [0.491252516,-0.132769930, 0.087742825],
            [0.539338687, 0.030052764, 0.117432841]
        ]
    }

    def to_dict(data):
        a = np.array(data)
        return {"x": a[:,0], "y": a[:,1], "yerr": a[:,2]}

    on  = {k: to_dict(v) for k,v in fa22_on.items()}
    off = {k: to_dict(v) for k,v in fa22_off.items()}

    fig, axes = plt.subplots(1, 3, figsize=(15,5))
    plt.suptitle("Fa22: Fiducial Cuts ON vs OFF", fontsize=16, y=0.97)
    fs = 13

    # ALU sinphi
    ax = axes[0]
    ax.errorbar(on["ALUsinphi"]["x"], on["ALUsinphi"]["y"],   yerr=on["ALUsinphi"]["yerr"],
                fmt="o", color=cmp_colors["on"], capsize=3, label="cuts ON")
    ax.errorbar(off["ALUsinphi"]["x"], off["ALUsinphi"]["y"], yerr=off["ALUsinphi"]["yerr"],
                fmt="o", color=cmp_colors["off"], capsize=3, label="cuts OFF")
    ax.set(xlim=(0,0.7), ylim=(-0.02,0.02),
           xlabel=r"$x_B$", ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    ax.axhline(0, linestyle="--", color="black");  ax.grid(True,linestyle="--",alpha=0.5)
    ax.legend(title="Config", frameon=True)

    # AUL sinphi and sin2phi
    ax = axes[1]
    # ON (circles)
    ax.errorbar(on["AULsinphi"]["x"], on["AULsinphi"]["y"], yerr=on["AULsinphi"]["yerr"],
                fmt="o", mfc="none", mec=cmp_colors["on"], capsize=3)
    ax.errorbar(on["AULsin2phi"]["x"], on["AULsin2phi"]["y"], yerr=on["AULsin2phi"]["yerr"],
                fmt="o", color=cmp_colors["on"], capsize=3)
    # OFF (squares)
    ax.errorbar(off["AULsinphi"]["x"], off["AULsinphi"]["y"], yerr=off["AULsinphi"]["yerr"],
                fmt="s", mfc="none", mec=cmp_colors["off"], capsize=3)
    ax.errorbar(off["AULsin2phi"]["x"],off["AULsin2phi"]["y"],yerr=off["AULsin2phi"]["yerr"],
                fmt="s", color=cmp_colors["off"], capsize=3)
    ax.set(xlim=(0,0.7), ylim=(-0.2,0.6),
           xlabel=r"$x_B$", ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    ax.axhline(0, linestyle="--", color="black");  ax.grid(True,linestyle="--",alpha=0.5)
    harm = ax.legend(
        handles=[
            Line2D([0],[0], marker='o', mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker='o', color='black', linestyle='', label='n=2')
        ],
        title="Harmonic", loc='upper right', frameon=True
    )
    ax.add_artist(harm)
    ax.legend(
        handles=[
            Line2D([0],[0], marker='o', color=cmp_colors["on"], linestyle='', label='cuts ON'),
            Line2D([0],[0], marker='s', color=cmp_colors["off"], linestyle='', label='cuts OFF')
        ],
        title="Config", loc='lower right', frameon=True
    )

    # ALL n=0 and cosphi
    ax = axes[2]
    ax.errorbar(on["ALL_n0"]["x"], on["ALL_n0"]["y"],   yerr=on["ALL_n0"]["yerr"],
                fmt="o", mfc="none", mec=cmp_colors["on"], capsize=3)
    ax.errorbar(on["ALLcosphi"]["x"], on["ALLcosphi"]["y"], yerr=on["ALLcosphi"]["yerr"],
                fmt="o", color=cmp_colors["on"], capsize=3)
    ax.errorbar(off["ALL_n0"]["x"], off["ALL_n0"]["y"], yerr=off["ALL_n0"]["yerr"],
                fmt="s", mfc="none", mec=cmp_colors["off"], capsize=3)
    ax.errorbar(off["ALLcosphi"]["x"],off["ALLcosphi"]["y"],yerr=off["ALLcosphi"]["yerr"],
                fmt="s", color=cmp_colors["off"], capsize=3)
    ax.set(xlim=(0,0.7), ylim=(-0.02,0.02),
           xlabel=r"$x_B$", ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    ax.axhline(0, linestyle="--", color="black");  ax.grid(True,linestyle="--",alpha=0.5)
    harm3 = ax.legend(
        handles=[
            Line2D([0],[0], marker='o', mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker='o', color='black', linestyle='', label='n=1')
        ],
        title="Harmonic", loc='upper right', frameon=True
    )
    ax.add_artist(harm3)
    ax.legend(
        handles=[
            Line2D([0],[0], marker='o', color=cmp_colors["on"], linestyle='', label='cuts ON'),
            Line2D([0],[0], marker='s', color=cmp_colors["off"], linestyle='', label='cuts OFF')
        ],
        title="Config", loc='lower right', frameon=True
    )

    plt.tight_layout(rect=[0,0,1,0.94])

    # Save but do not show
    cmp_file = os.path.join(out_dir, "fa22_cuts_on_vs_off.pdf")
    plt.savefig(cmp_file)
    print(f"Saved comparison plot to '{cmp_file}'")

# call it at script end
compare_fa22_versions()