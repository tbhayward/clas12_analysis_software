#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt

# Ensure output directory exists
os.makedirs("output/enpi+", exist_ok=True)

# ===================== Data (WITH A_tg turned ON) =====================

enpi_AULsinphi_on = [
    (-1.199181111, 0.044152285, 0.008751358), (-1.098543449, 0.030433354, 0.008364471),
    (-0.998883524, 0.031709080, 0.008408296), (-0.898912823, 0.034462384, 0.006377179),
    (-0.798568755, 0.026566700, 0.008667812), (-0.698809864, 0.027595580, 0.008566055),
    (-0.598626370, 0.012839146, 0.007449595), (-0.498526839, 0.021628335, 0.006477801),
    (-0.398650913, 0.029360333, 0.005293234), (-0.299064647, 0.025819033, 0.004363945),
    (-0.200361195, 0.022595616, 0.004226807), (-0.110121584, 0.016788042, 0.004935793),
]
enpi_AULsin2phi_on = [
    (-1.199181111, -0.108552020, 0.020379230), (-1.098543449, -0.083651810, 0.021115521),
    (-0.998883524, -0.135147222, 0.018408224), (-0.898912823, -0.079039866, 0.015652445),
    (-0.798568755, -0.109159950, 0.021737875), (-0.698809864, -0.111851880, 0.023115879),
    (-0.598626370, -0.057389535, 0.017827886), (-0.498526839, -0.077718447, 0.014222532),
    (-0.398650913, -0.044912422, 0.010774425), (-0.299064647, -0.017387861, 0.008385197),
    (-0.200361195, 0.004837571, 0.009423116),  (-0.110121584, -0.020247902, 0.013739176),
]
enpi_Atg_on = [
    (-1.199181111, -0.012793452, 0.023591008), (-1.098543449, 0.004175620, 0.024507994),
    (-0.998883524, 0.040753071, 0.021605454),  (-0.898912823, 0.018866257, 0.019857263),
    (-0.798568755, 0.037547198, 0.026235586),  (-0.698809864, 0.048891019, 0.026535179),
    (-0.598626370, -0.002182214, 0.021463908), (-0.498526839, 0.053014756, 0.015810257),
    (-0.398650913, 0.020564072, 0.011674460),  (-0.299064647, 0.006398819, 0.008616185),
    (-0.200361195, -0.006443835, 0.009862269), (-0.110121584, 0.010028340, 0.015131640),
]

lowxB_AULsinphi_on = [
    (-1.144968841, 0.053803942, 0.013434339), (-0.944704682, -0.022226282, 0.012094558),
    (-0.743862170, -0.046665066, 0.011426470), (-0.542210598, -0.052515848, 0.010767124),
    (-0.338175138, -0.045989046, 0.008144064), (-0.138113989, 0.019109932, 0.005530823),
]
lowxB_AULsin2phi_on = [
    (-1.144968841, -0.092980586, 0.043819335), (-0.944704682, -0.093403619, 0.025375769),
    (-0.743862170, -0.145105220, 0.025741785), (-0.542210598, -0.127116699, 0.034325640),
    (-0.338175138, -0.086784991, 0.018465639), (-0.138113989, -0.012565180, 0.014073559),
]
lowxB_Atg_on = [
    (-1.144968841, 0.049912118, 0.028503761), (-0.944704682, -0.003851290, 0.028675714),
    (-0.743862170, -0.031267062, 0.022155278), (-0.542210598, 0.058264292, 0.015807054),
    (-0.338175138, 0.027622693, 0.026227972), (-0.138113989, 0.020406214, 0.017383911),
]

highxB_AULsinphi_on = [
    (-1.145494332, 0.013279804, 0.012413866), (-0.945736523, 0.044797534, 0.008498857),
    (-0.745991387, 0.046421204, 0.008271592), (-0.551158219, 0.017595761, 0.007074191),
    (-0.381843164, 0.024728007, 0.009643565), (-0.240325241, 0.199998870, 0.011347087),
]
highxB_AULsin2phi_on = [
    (-1.145494332, -0.093085929, 0.018903874), (-0.945736523, -0.058980380, 0.019196957),
    (-0.745991387, -0.075057350, 0.017106973), (-0.551158219, -0.022290776, 0.014376402),
    (-0.381843164, 0.013443519, 0.017870449), (-0.240325241, -0.542910679, 0.074052279),
]
highxB_Atg_on = [
    (-1.145494332, 0.031123255, 0.017701595), (-0.945736523, 0.003716427, 0.016632532),
    (-0.745991387, 0.022760721, 0.013974358), (-0.551158219, -0.007066025, 0.012143044),
    (-0.381843164, -0.020106357, 0.019314443), (-0.240325241, -0.499999999, 0.002544515),
]

# ===================== Data (WITH A_tg turned OFF) =====================

enpi_AULsinphi_off = [
    (-1.199181111, 0.044093752, 0.009079273), (-1.098543449, 0.030053694, 0.007999242),
    (-0.998883524, 0.026308230, 0.007788450), (-0.898912823, 0.033068246, 0.006384161),
    (-0.798568755, 0.017585942, 0.006731740), (-0.698809864, 0.015322775, 0.005273293),
    (-0.598626370, 0.013430907, 0.004645720), (-0.498526839, 0.004885972, 0.004319493),
    (-0.398650913, 0.022804902, 0.003901649), (-0.299064647, 0.023768626, 0.003416027),
    (-0.200361195, 0.024272315, 0.003333052), (-0.110121584, 0.015303469, 0.004409024),
]
enpi_AULsin2phi_off = [
    (-1.199181111, -0.111799982, 0.020399582), (-1.098543449, -0.081553354, 0.017269733),
    (-0.998883524, -0.113543242, 0.014499052), (-0.898912823, -0.069388369, 0.011872551),
    (-0.798568755, -0.083272127, 0.010995075), (-0.698809864, -0.076260921, 0.010627099),
    (-0.598626370, -0.058909285, 0.009888374), (-0.498526839, -0.050099910, 0.008326457),
    (-0.398650913, -0.034044228, 0.007873420), (-0.299064647, -0.015547870, 0.007551487),
    (-0.200361195, 0.001486147, 0.009065971),  (-0.110121584, -0.014523537, 0.010276728),
]

lowxB_AULsinphi_off = [
    (-1.144968841, 0.037750296, 0.018800817), (-0.944704682, -0.022227775, 0.012090410),
    (-0.743862170, -0.046383860, 0.011760372), (-0.542210598, -0.067249039, 0.012226644),
    (-0.338175138, -0.048857560, 0.007713168), (-0.138113989, 0.015218672, 0.004489387),
]
lowxB_AULsin2phi_off = [
    (-1.144968841, -0.014033889, 0.048626163), (-0.944704682, -0.094860711, 0.022950509),
    (-0.743862170, -0.150330600, 0.027796755), (-0.542210598, -0.132423455, 0.021927933),
    (-0.338175138, -0.074196438, 0.012563805), (-0.138113989, -0.001110331, 0.009046579),
]

highxB_AULsinphi_off = [
    (-1.145494332, 0.028281724, 0.009035665), (-0.945736523, 0.045364777, 0.008123986),
    (-0.745991387, 0.040694119, 0.007691909), (-0.551158219, 0.019276513, 0.006434062),
    (-0.381843164, 0.028740387, 0.008815597), (-0.240325241, 0.200000000, 0.002776397),
]
highxB_AULsin2phi_off = [
    (-1.145494332, -0.084107471, 0.017544984), (-0.945736523, -0.059748406, 0.019243512),
    (-0.745991387, -0.061860181, 0.015150256), (-0.551158219, -0.023410365, 0.015467441),
    (-0.381843164, 0.015693699, 0.019334087),  (-0.240325241, -0.542981435, 0.036034762),
]

# ===================== Helpers =====================

def to_xyyerr(triplets):
    """Convert list of (x_raw, y, err) -> arrays with x = -x_raw (plot versus -t)."""
    arr = np.array(triplets, dtype=float)
    x = -arr[:, 0]  # -t on horizontal axis
    y = arr[:, 1]
    e = arr[:, 2]
    idx = np.argsort(x)
    return x[idx], y[idx], e[idx]

def add_zero_line(ax):
    ax.axhline(0.0, lw=1.0, alpha=0.5)

# ===================== Plotting =====================

fig, axes = plt.subplots(3, 3, figsize=(14, 10), sharex=False, sharey=False)
plt.subplots_adjust(wspace=0.25, hspace=0.28)

# Column titles
col_titles = [r"$A_{UL}^{\sin\phi}$", r"$A_{UL}^{\sin 2\phi}$", r"$A_{tg}$"]
for j, title in enumerate(col_titles):
    axes[0, j].set_title(title, fontsize=13, pad=8)

# Row labels
row_labels = [r"$x_B\in[0.10,\,0.60]$ (integrated)",
              r"$x_B\in[0.10,\,0.25]$ (low)",
              r"$x_B\in[0.45,\,0.60]$ (high)"]

# ------------- Row 1: integrated over x_B -------------
# Col 1: AUL^sinφ
ax = axes[0, 0]
x_on, y_on, e_on = to_xyyerr(enpi_AULsinphi_on)
x_off, y_off, e_off = to_xyyerr(enpi_AULsinphi_off)
ax.errorbar(x_on,  y_on,  yerr=e_on,  fmt='o', linestyle='none', label='with $A_{tg}$', capsize=3)
ax.errorbar(x_off, y_off, yerr=e_off, fmt='s', linestyle='none', label='no $A_{tg}$', capsize=3)
add_zero_line(ax)
ax.set_ylabel("Amplitude")
ax.legend(frameon=False, fontsize=9)
ax.text(0.02, 0.95, row_labels[0], transform=ax.transAxes, va='top', fontsize=11)

# Col 2: AUL^sin2φ
ax = axes[0, 1]
x_on, y_on, e_on = to_xyyerr(enpi_AULsin2phi_on)
x_off, y_off, e_off = to_xyyerr(enpi_AULsin2phi_off)
ax.errorbar(x_on,  y_on,  yerr=e_on,  fmt='o', linestyle='none', label='with $A_{tg}$', capsize=3)
ax.errorbar(x_off, y_off, yerr=e_off, fmt='s', linestyle='none', label='no $A_{tg}$', capsize=3)
add_zero_line(ax)
ax.legend(frameon=False, fontsize=9)

# Col 3: A_tg
ax = axes[0, 2]
x_atg, y_atg, e_atg = to_xyyerr(enpi_Atg_on)
ax.errorbar(x_atg, y_atg, yerr=e_atg, fmt='o', linestyle='none', label='$A_{tg}$ fit', capsize=3)
add_zero_line(ax)
ax.legend(frameon=False, fontsize=9)

# ------------- Row 2: low x_B -------------
# Col 1: AUL^sinφ
ax = axes[1, 0]
x_on, y_on, e_on = to_xyyerr(lowxB_AULsinphi_on)
x_off, y_off, e_off = to_xyyerr(lowxB_AULsinphi_off)
ax.errorbar(x_on,  y_on,  yerr=e_on,  fmt='o', linestyle='none', label='with $A_{tg}$', capsize=3)
ax.errorbar(x_off, y_off, yerr=e_off, fmt='s', linestyle='none', label='no $A_{tg}$', capsize=3)
add_zero_line(ax)
ax.set_ylabel("Amplitude")
ax.legend(frameon=False, fontsize=9)
ax.text(0.02, 0.95, row_labels[1], transform=ax.transAxes, va='top', fontsize=11)

# Col 2: AUL^sin2φ
ax = axes[1, 1]
x_on, y_on, e_on = to_xyyerr(lowxB_AULsin2phi_on)
x_off, y_off, e_off = to_xyyerr(lowxB_AULsin2phi_off)
ax.errorbar(x_on,  y_on,  yerr=e_on,  fmt='o', linestyle='none', label='with $A_{tg}$', capsize=3)
ax.errorbar(x_off, y_off, yerr=e_off, fmt='s', linestyle='none', label='no $A_{tg}$', capsize=3)
add_zero_line(ax)
ax.legend(frameon=False, fontsize=9)

# Col 3: A_tg
ax = axes[1, 2]
x_atg, y_atg, e_atg = to_xyyerr(lowxB_Atg_on)
ax.errorbar(x_atg, y_atg, yerr=e_atg, fmt='o', linestyle='none', label='$A_{tg}$ fit', capsize=3)
add_zero_line(ax)
ax.legend(frameon=False, fontsize=9)

# ------------- Row 3: high x_B -------------
# Col 1: AUL^sinφ
ax = axes[2, 0]
x_on, y_on, e_on = to_xyyerr(highxB_AULsinphi_on)
x_off, y_off, e_off = to_xyyerr(highxB_AULsinphi_off)
ax.errorbar(x_on,  y_on,  yerr=e_on,  fmt='o', linestyle='none', label='with $A_{tg}$', capsize=3)
ax.errorbar(x_off, y_off, yerr=e_off, fmt='s', linestyle='none', label='no $A_{tg}$', capsize=3)
add_zero_line(ax)
ax.set_ylabel("Amplitude")
ax.set_xlabel(r"$-t\;(\mathrm{GeV}^2)$")
ax.legend(frameon=False, fontsize=9)
ax.text(0.02, 0.95, row_labels[2], transform=ax.transAxes, va='top', fontsize=11)

# Col 2: AUL^sin2φ
ax = axes[2, 1]
x_on, y_on, e_on = to_xyyerr(highxB_AULsin2phi_on)
x_off, y_off, e_off = to_xyyerr(highxB_AULsin2phi_off)
ax.errorbar(x_on,  y_on,  yerr=e_on,  fmt='o', linestyle='none', label='with $A_{tg}$', capsize=3)
ax.errorbar(x_off, y_off, yerr=e_off, fmt='s', linestyle='none', label='no $A_{tg}$', capsize=3)
add_zero_line(ax)
ax.set_xlabel(r"$-t\;(\mathrm{GeV}^2)$")
ax.legend(frameon=False, fontsize=9)

# Col 3: A_tg
ax = axes[2, 2]
x_atg, y_atg, e_atg = to_xyyerr(highxB_Atg_on)
ax.errorbar(x_atg, y_atg, yerr=e_atg, fmt='o', linestyle='none', label='$A_{tg}$ fit', capsize=3)
add_zero_line(ax)
ax.set_xlabel(r"$-t\;(\mathrm{GeV}^2)$")
ax.legend(frameon=False, fontsize=9)

# Tidy axes: grid, ticks, fixed y-limits, reasonable x-limits
for i in range(3):
    for j in range(3):
        axes[i, j].grid(alpha=0.25)
        axes[i, j].tick_params(direction='in')
        axes[i, j].set_ylim(-0.2, 0.2)  # standardized y-axis across all panels
        if i < 2:
            axes[i, j].set_xlabel("")
        # Auto x-limits from plotted data
        xs = []
        for line in axes[i, j].lines:
            xd = line.get_xdata()
            if len(xd): xs.extend(list(xd))
        if xs:
            xmin, xmax = min(xs), max(xs)
            pad = 0.03 * (xmax - xmin if xmax > xmin else 1.0)
            axes[i, j].set_xlim(xmin - pad, xmax + pad)

fig.suptitle(r"UL amplitudes vs $-t$ (with and without tangential leakage term)", fontsize=14, y=0.99)
plt.savefig("output/enpi+/GE_UL_3x3_panels.pdf", dpi=300, bbox_inches="tight")
plt.show()