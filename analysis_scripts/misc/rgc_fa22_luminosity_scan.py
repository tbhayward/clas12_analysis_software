#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

# Data: beam current (nA), counts, charge (nC)
currents = np.array([2, 4, 6, 8, 10])
counts   = np.array([656206, 1321205, 1729427, 1167246, 1032345])
charges  = np.array([5969.59, 12153.266, 16195.799, 11276.1231, 9217.6748])

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

# Make sure output directory exists
os.makedirs("output", exist_ok=True)

# Plot
plt.errorbar(currents, counts_per_nC, yerr=uncertainties, fmt='o', label='Data')
x_fit = np.linspace(currents.min(), currents.max(), 100)
y_fit = m * x_fit + b
plt.plot(x_fit, y_fit, '-', label=f'Fit: y = {m:.2f} x + {b:.2f}')

plt.xlabel("nA")
plt.ylabel("counts/nC")
plt.legend(loc="best")

plt.tight_layout()
plt.savefig("output/rgc_fa22_luminosity_scan.pdf")
plt.close()