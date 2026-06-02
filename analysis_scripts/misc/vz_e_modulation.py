#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt


def modulation_label(modulation):
    labels = {
        "ALUsinphi": r"$F_{LU}^{\sin\phi}/F_{UU}$",
        "AULsinphi": r"$F_{UL}^{\sin\phi}/F_{UU}$",
        "AULsin2phi": r"$F_{UL}^{\sin 2\phi}/F_{UU}$",
        "ALL": r"$F_{LL}/F_{UU}$",
        "ALLcosphi": r"$F_{LL}^{\cos\phi}/F_{UU}$",
    }

    return labels[modulation]


def make_vz_e_modulation_plot():
    output_path = "output/vz_e_modulations.png"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    data = {
        "ALUsinphi": np.array([
            [-4.512208106, 0.103969702, 0.008360384],
            [-2.789373615, 0.107822857, 0.006620113],
            [-0.897536588, 0.097091414, 0.006873885],
        ], dtype=float),

        "AULsinphi": np.array([
            [-4.512208106, 0.078894027, 0.007775402],
            [-2.789373615, 0.068056361, 0.006029922],
            [-0.897536588, 0.050513240, 0.006252694],
        ], dtype=float),

        "AULsin2phi": np.array([
            [-4.512208106, -0.131050328, 0.018504594],
            [-2.789373615, -0.130111872, 0.012969811],
            [-0.897536588, -0.151914943, 0.013539135],
        ], dtype=float),

        "ALL": np.array([
            [-4.512208106, 0.612054281, 0.017540018],
            [-2.789373615, 0.620648227, 0.013598041],
            [-0.897536588, 0.606664748, 0.014639742],
        ], dtype=float),

        "ALLcosphi": np.array([
            [-4.512208106, 0.097070145, 0.079934250],
            [-2.789373615, 0.115280769, 0.063959175],
            [-0.897536588, 0.208809058, 0.056005670],
        ], dtype=float),
    }

    modulations = [
        "ALUsinphi",
        "AULsinphi",
        "AULsin2phi",
        "ALL",
        "ALLcosphi",
    ]

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 8.0),
        sharex=True,
    )

    axes_flat = axes.flatten()

    for i_mod, modulation in enumerate(modulations):
        ax = axes_flat[i_mod]
        arr = data[modulation]

        vz_e = arr[:, 0]
        asymmetry = arr[:, 1]
        asymmetry_err = arr[:, 2]

        ax.errorbar(
            vz_e,
            asymmetry,
            yerr=asymmetry_err,
            fmt="o",
            markersize=5,
            capsize=3,
            linewidth=1,
        )

        ax.axhline(0.0, linestyle="--", linewidth=1)
        ax.grid(True, alpha=0.3)

        ax.set_title(modulation_label(modulation), fontsize=15)
        ax.set_xlabel(r"$v_{z}^{e}$ (cm)", fontsize=13)
        ax.set_ylabel("asymmetry", fontsize=13)

        y_min = np.min(asymmetry - asymmetry_err)
        y_max = np.max(asymmetry + asymmetry_err)
        y_pad = 0.20 * (y_max - y_min)

        if y_pad <= 0.0:
            y_pad = 0.01
        #endif

        ax.set_ylim(y_min - y_pad, y_max + y_pad)
    #endfor

    axes_flat[5].axis("off")

    fig.suptitle(
        r"Exclusive $e p \rightarrow e n \pi^{+}$ polarized modulations vs. $v_{z}^{e}$",
        fontsize=17,
    )

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Saved plot to {output_path}")


def main():
    make_vz_e_modulation_plot()


if __name__ == "__main__":
    main()
#endif