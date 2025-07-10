#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Ensure output directory exists
os.makedirs('output', exist_ok=True)

# Data sets
xchi2FitsALUsinphi = [
    [
        (0.139608167, 0.036519444, 0.004440912),
        (0.223826090, 0.047197386, 0.004232517),
        (0.323794485, 0.048494028, 0.007532405),
        (0.423461067, 0.047274189, 0.008662584),
        (0.518591309, 0.060445966, 0.014232318),
        (0.607365519, 0.036538084, 0.038878480),
    ],
    [
        (0.139608167, 0.036748334, 0.004617308),
        (0.223826090, 0.048331934, 0.004096731),
        (0.323794485, 0.054207700, 0.006117512),
        (0.423461067, 0.045095539, 0.009356716),
        (0.518591309, 0.054181510, 0.014506165),
        (0.607365519, 0.031391103, 0.037215555),
    ],
    [
        (0.139608167, 0.031868239, 0.005041277),
        (0.223826090, 0.047092432, 0.004009034),
        (0.323794485, 0.044737175, 0.008312746),
        (0.423461067, 0.047302695, 0.008446926),
        (0.518591309, 0.063933250, 0.013918055),
        (0.607365519, 0.009842388, 0.038626854),
    ],
]

# LaTeX labels
labels = [r'$ep\pi^{0}$', r'$\gamma_{1}$', r'$\gamma_{2}$']

# Plotting loop
for data, label in zip(xchi2FitsALUsinphi, labels):
    xs = [point[0] for point in data]
    ys = [point[1] for point in data]
    yerrs = [point[2] for point in data]
    plt.errorbar(xs, ys, yerr=yerrs, fmt='o', label=label)
#endfor

# Axis labels and limits
plt.xlabel(r'$x_{B}$')
plt.ylabel(r'$A_{LU}$')
plt.xlim(0, 0.7)
plt.ylim(-0.1, 0.1)

# Legend and layout
plt.legend()
plt.tight_layout()

# Save figure
plt.savefig('output/eppi0_ALU_gamma_study.pdf')
plt.close()