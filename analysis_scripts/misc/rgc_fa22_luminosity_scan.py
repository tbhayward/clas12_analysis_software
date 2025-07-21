#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

# Data: beam current (nA), counts, charge (nC)
# Removed the 10 nA point
currents = np.array([2, 4, 6, 8])
counts   = np.array([656206, 1321205, 1729427, 1167246])
charges  = np.array([5969.59, 12153.266, 16195.799, 11276.1231])

# Compute counts per nC and its uncertainty (sqrt(counts)/charge)
counts_per_nC = []
uncertainties = []
for cnt, q in zip(counts, charges):
    counts_per_nC.append(cnt / q)
    uncertainties.append(np.sqrt(cnt) / q)
#endfor

counts_per_nC = np.array(counts_per_nC)
uncertainties = np.array(uncertainties)

# Fit a straight line y = m*x + b
m, b = np.polyfit(currents, counts_per_nC, 1)

# Compute fitted value at 0 nA and percentage drops
y0 = b  # value at 0 nA
drop_percentages = []
for I in currents:
    yI = m * I + b
    drop_pct = (y0 - yI) / y0 * 100
    drop_percentages.append(drop_pct)
#endfor

# Prepare table data: first row is fitted value, then drops
row_labels = [
    "0 nA (counts/nC)",
    "2 nA drop (%)",
    "4 nA drop (%)",
    "6 nA drop (%)",
    "8 nA drop (%)",
]
table_data = [[f"{y0:.2f}"]] + [[f"{dp:.2f}%"] for dp in drop_percentages]

# Ensure output directory exists
os.makedirs("output", exist_ok=True)

# Create plot
fig, ax = plt.subplots()
ax.errorbar(currents, counts_per_nC, yerr=uncertainties, fmt='o', label='Data')
x_fit = np.linspace(currents.min(), currents.max(), 100)
y_fit = m * x_fit + b
ax.plot(x_fit, y_fit, '-', label=f'Fit: y = {m:.2f} x + {b:.2f}')

# Add table in the top right
table = ax.table(
    cellText=table_data,
    rowLabels=row_labels,
    colLabels=["Value"],
    cellLoc='center',
    rowLoc='center',
    loc='upper right',
)
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.5)

# Labels and legend
ax.set_xlabel("nA")
ax.set_ylabel("counts/nC")
ax.legend(loc="lower left")

plt.tight_layout()
plt.savefig("output/rgc_fa22_luminosity_scan.pdf")
plt.close()