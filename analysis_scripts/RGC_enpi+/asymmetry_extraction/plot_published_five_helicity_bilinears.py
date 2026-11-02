#!/usr/bin/env python3
"""Plot five helicity bilinears derived from the hard-coded analysis-note tables.

The input values are Tables 25, 28, 31, and 34 of the supplied RGC analysis
note.  Each entry is (central value, statistical error, point-to-point
systematic error).  The note does not publish the five-observable covariance,
so uncertainties of sums/differences are propagated diagonally and labelled as
such.  Correlated scale uncertainties are not included, matching the tables.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from helicity_bilinear_transforms import (
    PUBLISHED_FIVE_LABELS,
    PUBLISHED_FIVE_PARAMETERS,
    published_five_diagonal_uncertainties,
    published_five_from_ratios,
)


XB_GROUPS = (
    (0.10, 0.25, "#0072B2", "o"),
    (0.25, 0.35, "#D55E00", "s"),
    (0.35, 0.45, "#009E73", "^"),
    (0.45, 0.60, "#CC79A7", "D"),
)

# Columns inside each row:
# -t', then for (lu1, ul1, ul2, ll0, ll1), repeating value, stat, ptp syst.
NOTE_TABLES = (
    (
        (0.129, 0.130, .011, .011, -0.031, .012, .033, -0.043, .023, .037, 0.333, .029, .031, 0.419, .041, .031),
        (0.340, 0.119, .016, .017, -0.194, .021, .041, -0.222, .035, .040, 0.395, .042, .024, 0.401, .055, .024),
        (0.544, 0.069, .020, .024, -0.184, .027, .023, -0.281, .043, .055, 0.463, .067, .020, 0.400, .098, .062),
        (0.745, 0.064, .018, .026, -0.078, .023, .013, -0.270, .050, .034, 0.345, .061, .073, -0.020, .123, .209),
        (0.945, 0.067, .028, .063,  0.026, .038, .011, -0.276, .066, .088, 0.393, .064, .094, -0.029, .176, .132),
        (1.145, 0.008, .027, .077,  0.113, .046, .076, -0.120, .120, .156, 0.369, .092, .049, 0.162, .160, .343),
    ),
    (
        (0.134, 0.151, .010, .006,  0.039, .008, .016, -0.097, .016, .013, 0.548, .022, .026, 0.247, .062, .018),
        (0.341, 0.138, .017, .013, -0.021, .012, .036, -0.162, .022, .035, 0.652, .031, .045, 0.485, .065, .171),
        (0.543, 0.107, .017, .017,  0.003, .016, .024, -0.240, .030, .071, 0.691, .041, .032, 0.227, .116, .150),
        (0.744, 0.120, .020, .016,  0.006, .018, .026, -0.255, .038, .041, 0.576, .051, .079, 0.152, .107, .081),
        (0.946, 0.095, .025, .018,  0.100, .025, .035, -0.222, .047, .127, 0.610, .057, .056, -0.039, .191, .320),
        (1.145, 0.143, .025, .053,  0.151, .028, .039, -0.213, .049, .034, 0.595, .061, .089, -0.074, .156, .103),
    ),
    (
        (0.141, 0.147, .011, .028, 0.078, .008, .006, -0.072, .017, .009, 0.692, .022, .030, 0.114, .088, .051),
        (0.343, 0.167, .014, .009, 0.058, .011, .012, -0.133, .022, .028, 0.747, .030, .050, 0.201, .097, .072),
        (0.544, 0.134, .016, .022, 0.043, .013, .012, -0.180, .025, .024, 0.719, .035, .030, 0.146, .108, .027),
        (0.745, 0.158, .020, .024, 0.090, .015, .015, -0.185, .028, .029, 0.696, .041, .036, -0.130, .118, .142),
        (0.944, 0.154, .023, .011, 0.112, .019, .018, -0.262, .042, .029, 0.738, .050, .053, 0.029, .140, .090),
        (1.145, 0.111, .026, .043, 0.111, .023, .031, -0.315, .048, .051, 0.728, .065, .053, -0.061, .174, .129),
    ),
    (
        (0.149, 0.119, .015, .015, 0.061, .011, .015, -0.070, .025, .022, 0.838, .030, .029, 0.175, .169, .077),
        (0.345, 0.172, .017, .025, 0.102, .013, .029, -0.104, .028, .023, 0.876, .035, .058, -0.102, .164, .400),
        (0.546, 0.132, .019, .016, 0.071, .016, .021, -0.137, .030, .031, 0.879, .040, .027, -0.039, .233, .481),
        (0.746, 0.092, .022, .030, 0.077, .019, .025, -0.175, .033, .021, 0.885, .049, .050, -0.253, .238, .113),
        (0.946, 0.161, .025, .033, 0.119, .023, .024, -0.210, .041, .065, 0.806, .051, .061, -0.311, .248, .427),
        (1.147, 0.120, .026, .025, 0.085, .022, .027, -0.149, .053, .051, 0.820, .061, .066, 0.291, .193, .303),
    ),
)


def unpack_group(group: tuple[tuple[float, ...], ...]):
    data = np.asarray(group, dtype=float)
    x = data[:, 0]
    values = data[:, [1, 4, 7, 10, 13]]
    stat = data[:, [2, 5, 8, 11, 14]]
    syst = data[:, [3, 6, 9, 12, 15]]
    return (
        x,
        published_five_from_ratios(values),
        published_five_diagonal_uncertainties(stat),
        published_five_diagonal_uncertainties(syst),
    )


def draw_group(group_index: int, output_dir: Path) -> Path:
    x, values, stat, syst = unpack_group(NOTE_TABLES[group_index])
    low, high, color, marker = XB_GROUPS[group_index]
    fig, axes = plt.subplots(2, 3, figsize=(20, 12), sharex=True)
    axes = axes.ravel()
    for panel, name in enumerate(PUBLISHED_FIVE_PARAMETERS):
        axis = axes[panel]
        axis.axhline(0.0, color="0.55", ls="--", lw=1.4)
        axis.bar(x, 2.0 * syst[:, panel], bottom=-syst[:, panel], width=.045,
                 color=color, alpha=.16, edgecolor="none",
                 label="P-to-p syst. (diagonal propagation)")
        axis.errorbar(x, values[:, panel], yerr=stat[:, panel], fmt=marker,
                      ms=8, mfc="white", mec=color, mew=2.2, ecolor=color,
                      capsize=4, label="Derived RGC bilinear (stat.)")
        axis.set_title(PUBLISHED_FIVE_LABELS[name], fontsize=11, pad=9)
        axis.set_ylabel("normalized bilinear")
        axis.text(.04, .91, f"({chr(97 + panel)})", transform=axis.transAxes,
                  fontsize=15, fontweight="bold")
        axis.minorticks_on()
        axis.tick_params(direction="in", top=True, right=True, which="both")
    axes[-1].axis("off")
    for axis in axes[3:5]:
        axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^{2})$", fontsize=14)
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="center", frameon=False, fontsize=12)
    axes[-1].text(
        .5, .25,
        "Input: analysis-note published five\n"
        "UL/LL labels retain the lab-frame provenance\n"
        "No unpolarized cos(phi) or cos(2phi) input",
        ha="center", va="center", transform=axes[-1].transAxes, fontsize=12,
    )
    fig.suptitle(
        rf"CLAS12 $ep\rightarrow e(n)\pi^+$ normalized helicity bilinears: ${low:.2f}<x_B<{high:.2f}$",
        fontsize=21,
    )
    fig.text(
        .5, .015,
        r"$T_{\mu,\lambda}^{0\sigma}$: photon helicity $\mu$, initial-proton helicity $\lambda$, "
        r"pion helicity 0, recoil-neutron helicity $\sigma$; sum over $\sigma=\pm\frac{1}{2}$.",
        ha="center", fontsize=11,
    )
    fig.tight_layout(rect=(.03, .05, .98, .95))
    path = output_dir / f"published_five_helicity_bilinears_xB_{group_index + 1}.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def write_csv(output_dir: Path) -> Path:
    path = output_dir / "published_five_helicity_bilinears.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        header = ["xB_min", "xB_max", "mean_minus_tprime_gev2"]
        for name in PUBLISHED_FIVE_PARAMETERS:
            header.extend((name, f"{name}_stat_diagonal", f"{name}_ptp_syst_diagonal"))
        writer.writerow(header)
        for group_index, group in enumerate(NOTE_TABLES):
            x, values, stat, syst = unpack_group(group)
            low, high, _, _ = XB_GROUPS[group_index]
            for row in range(len(x)):
                output = [low, high, x[row]]
                for column in range(5):
                    output.extend((values[row, column], stat[row, column], syst[row, column]))
                writer.writerow(output)
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path,
                        default=Path("output/published_five_helicity_bilinears"))
    args = parser.parse_args()
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    csv_path = write_csv(output)
    plots = [draw_group(index, output) for index in range(4)]
    print(f"Wrote {csv_path}")
    for path in plots:
        print(f"Wrote {path}")
    print("WARNING: note-table uncertainty propagation is diagonal because the published covariance is unavailable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
