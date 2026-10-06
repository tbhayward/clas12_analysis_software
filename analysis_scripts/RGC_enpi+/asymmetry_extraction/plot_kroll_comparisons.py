#!/usr/bin/env python3
"""Publication-style comparison of final RGC e n pi+ ratios with Kroll predictions.

Run this *after* extract_structure_function_ratios.py.  By default the script reads

    output/nominal/tables/structure_function_ratios.csv

and writes PNG files to

    output/kroll_comparisons/

The Kroll predictions supplied in October 2026 are hard-coded below.  They are
virtual-photon-frame quantities.  The RGC UL/LL measurements remain in the lab-
frame convention used by the analysis; no frame conversion is applied here.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path("output/nominal/tables/structure_function_ratios.csv")
DEFAULT_OUTPUT_DIR = Path("output/kroll_comparisons")
DPI = 300

XB_BINS = (
    (0.10, 0.25),
    (0.25, 0.35),
    (0.35, 0.45),
    (0.45, 0.60),
)

XB_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
XB_MARKERS = ("o", "s", "^", "D")

OBSERVABLES = (
    "lu1",
    "ul1",
    "ul2",
    "ll0",
    "ll1",
)

OBSERVABLE_LABELS = {
    "lu1": r"$F_{LU}^{\sin\phi}/F_{UU}$",
    "ul1": r"$F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$",
    "ul2": r"$F_{UL,\mathrm{lab}}^{\sin2\phi}/F_{UU}$",
    "ll0": r"$F_{LL,\mathrm{lab}}/F_{UU}$",
    "ll1": r"$F_{LL,\mathrm{lab}}^{\cos\phi}/F_{UU}$",
}

# Kroll labels deliberately omit "lab": these are virtual-photon-frame predictions.
KROLL_LABELS = {
    "lu1": r"$F_{LU}^{\sin\phi}/F_{UU}$",
    "ul1": r"$F_{UL}^{\sin\phi}/F_{UU}$",
    "ul2": r"$F_{UL}^{\sin2\phi}/F_{UU}$",
    "ll0": r"$F_{LL}/F_{UU}$",
    "ll1": r"$F_{LL}^{\cos\phi}/F_{UU}$",
}

SINGLE_SPIN = frozenset(("lu1", "ul1", "ul2"))

# The table supplied by Peter Kroll: xB, -t', Q2, W, epsilon,
# F_LU^sin(phi)/F_UU, F_UL^sin(phi)/F_UU, F_UL^sin(2phi)/F_UU,
# F_LL/F_UU, F_LL^cos(phi)/F_UU.
KROLL_ROWS = (
    (0.2, 0.2, 2.13, 3.066, 0.719, 0.192, -0.230, -0.102, 0.929, -0.068),
    (0.2, 0.4, 2.13, 3.066, 0.712, 0.136, -0.160, -0.104, 0.933, -0.039),
    (0.2, 0.6, 2.13, 3.066, 0.700, 0.118, -0.133, -0.109, 0.942, -0.032),
    (0.2, 0.8, 2.13, 3.066, 0.695, 0.114, -0.123, -0.119, 0.954, -0.035),
    (0.2, 1.0, 2.13, 3.066, 0.698, 0.119, -0.122, -0.133, 0.963, -0.042),
    (0.2, 1.2, 2.13, 3.066, 0.694, 0.129, -0.129, -0.152, 0.969, -0.053),
    (0.3, 0.2, 2.52, 2.600, 0.823, 0.159, -0.203, -0.082, 0.939, -0.007),
    (0.3, 0.4, 2.52, 2.600, 0.825, 0.122, -0.155, -0.090, 0.938, -0.010),
    (0.3, 0.6, 2.52, 2.600, 0.828, 0.105, -0.131, -0.094, 0.942, -0.016),
    (0.3, 0.8, 2.52, 2.600, 0.822, 0.097, -0.120, -0.098, 0.948, -0.024),
    (0.3, 1.0, 2.52, 2.600, 0.819, 0.095, -0.116, -0.105, 0.954, -0.033),
    (0.3, 1.2, 2.52, 2.600, 0.819, 0.096, -0.115, -0.112, 0.959, -0.043),
    (0.4, 0.2, 3.01, 2.323, 0.857, 0.111, -0.154, -0.057, 0.949, -0.023),
    (0.4, 0.4, 3.01, 2.323, 0.856, 0.097, -0.133, -0.071, 0.945, -0.024),
    (0.4, 0.6, 3.01, 2.323, 0.858, 0.089, -0.118, -0.077, 0.948, -0.026),
    (0.4, 0.8, 3.01, 2.323, 0.857, 0.085, -0.109, -0.083, 0.954, -0.029),
    (0.4, 1.0, 3.01, 2.323, 0.853, 0.085, -0.105, -0.089, 0.960, -0.034),
    (0.4, 1.2, 3.01, 2.323, 0.849, 0.087, -0.103, -0.096, 0.966, -0.039),
    (0.5, 0.2, 4.24, 2.263, 0.829, 0.081, -0.117, -0.036, 0.961, -0.009),
    (0.5, 0.4, 4.24, 2.263, 0.825, 0.083, -0.117, -0.051, 0.956, -0.016),
    (0.5, 0.6, 4.24, 2.263, 0.824, 0.082, -0.111, -0.060, 0.955, -0.022),
    (0.5, 0.8, 4.24, 2.263, 0.824, 0.082, -0.107, -0.067, 0.957, -0.028),
    (0.5, 1.0, 4.24, 2.263, 0.822, 0.083, -0.105, -0.072, 0.960, -0.034),
    (0.5, 1.2, 4.24, 2.263, 0.821, 0.086, -0.104, -0.078, 0.964, -0.039),
)


def configure_style() -> None:
    """Use a clean, journal-friendly Matplotlib style without requiring LaTeX."""
    mpl.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": DPI,
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 11.5,
        "axes.labelsize": 13,
        "axes.titlesize": 12.5,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 10,
        "axes.linewidth": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 5.0,
        "ytick.major.size": 5.0,
        "xtick.minor.size": 2.8,
        "ytick.minor.size": 2.8,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "lines.linewidth": 1.8,
        "axes.unicode_minus": False,
    })


def kroll_frame() -> pd.DataFrame:
    columns = ("xB", "minus_tprime", "Q2", "W", "epsilon", *OBSERVABLES)
    return pd.DataFrame(KROLL_ROWS, columns=columns)


def validate_measurements(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"bin_number", "mean_xB", "mean_minus_tprime_gev2"}
    for observable in OBSERVABLES:
        required.add(observable)
        required.add(f"{observable}_stat")
        required.add(f"{observable}_point_to_point_systematic")
    # endfor

    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(
            "Input CSV is missing required final-result columns:\n  "
            + "\n  ".join(missing)
        )
    # endif

    output = frame.copy().sort_values("bin_number").reset_index(drop=True)
    if len(output) != 24:
        raise ValueError(f"Expected 24 analysis bins, found {len(output)} rows.")
    # endif

    output["x_index_plot"] = (output["bin_number"].astype(int) - 1) // 6
    return output


def add_systematic_boxes(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    sys: np.ndarray,
    color: str,
    width: float = 0.055,
    alpha: float = 0.20,
    zorder: float = 2.0,
) -> None:
    """Draw point-to-point systematic uncertainties as centered boxes."""
    for xi, yi, si in zip(x, y, sys):
        if not (np.isfinite(xi) and np.isfinite(yi) and np.isfinite(si)):
            continue
        # endif
        rect = Rectangle(
            (xi - 0.5 * width, yi - si),
            width,
            2.0 * si,
            facecolor=color,
            edgecolor="none",
            alpha=alpha,
            zorder=zorder,
        )
        ax.add_patch(rect)
    # endfor


def draw_measurements(
    ax: plt.Axes,
    subset: pd.DataFrame,
    observable: str,
    color: str,
    marker: str,
    label: str | None = None,
) -> None:
    x = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
    y = subset[observable].to_numpy(dtype=float)
    stat = subset[f"{observable}_stat"].to_numpy(dtype=float)
    sys = subset[f"{observable}_point_to_point_systematic"].to_numpy(dtype=float)

    add_systematic_boxes(ax, x, y, sys, color=color)
    ax.errorbar(
        x,
        y,
        yerr=stat,
        fmt=marker,
        ms=6.2,
        mfc="white",
        mec=color,
        mew=1.35,
        ecolor=color,
        elinewidth=1.25,
        capsize=2.4,
        capthick=1.1,
        linestyle="none",
        color=color,
        label=label,
        zorder=4,
    )


def draw_kroll_curve(
    ax: plt.Axes,
    predictions: pd.DataFrame,
    observable: str,
    color: str,
    label: str | None = None,
) -> None:
    x = predictions["minus_tprime"].to_numpy(dtype=float)
    y = predictions[observable].to_numpy(dtype=float)

    # At the forward limit the sine modulations vanish.  Add that physical
    # endpoint only for the three single-spin sine observables, as requested.
    if observable in SINGLE_SPIN:
        x = np.concatenate(([0.0], x))
        y = np.concatenate(([0.0], y))
    # endif

    ax.plot(
        x,
        y,
        color=color,
        lw=2.0,
        linestyle="-",
        label=label,
        zorder=3,
    )


def format_axis(ax: plt.Axes, observable: str, show_xlabel: bool = True) -> None:
    ax.axhline(0.0, color="0.55", lw=0.8, linestyle="--", zorder=0)
    ax.set_xlim(-0.035, 1.285)
    ax.set_xticks(np.arange(0.0, 1.21, 0.2))
    ax.minorticks_on()
    ax.set_ylabel(OBSERVABLE_LABELS[observable])
    if show_xlabel:
        ax.set_xlabel(r"$-t'$ (GeV$^2$)")
    # endif
    ax.grid(False)


def robust_y_limits(
    measurement_subsets: list[pd.DataFrame],
    prediction_subsets: list[pd.DataFrame],
    observable: str,
) -> tuple[float, float]:
    values: list[float] = [0.0]
    for subset in measurement_subsets:
        y = subset[observable].to_numpy(dtype=float)
        stat = subset[f"{observable}_stat"].to_numpy(dtype=float)
        sys = subset[f"{observable}_point_to_point_systematic"].to_numpy(dtype=float)
        total_extent = stat + sys
        values.extend((y - total_extent)[np.isfinite(y - total_extent)].tolist())
        values.extend((y + total_extent)[np.isfinite(y + total_extent)].tolist())
    # endfor
    for subset in prediction_subsets:
        y = subset[observable].to_numpy(dtype=float)
        values.extend(y[np.isfinite(y)].tolist())
    # endfor

    lo = float(np.min(values))
    hi = float(np.max(values))
    span = max(hi - lo, 0.12)
    pad = 0.12 * span
    lo -= pad
    hi += pad

    # Keep the near-unity LL constant term visually resolved instead of
    # needlessly forcing zero into that panel.
    if observable == "ll0" and lo > 0.25:
        lo = max(0.0, lo)
    # endif
    return lo, hi


def info_panel(ax: plt.Axes, x_index: int | None = None) -> None:
    ax.axis("off")
    handles = [
        Line2D([], [], marker="o", linestyle="none", mfc="white", mec="black",
               mew=1.3, ms=6.5, label="RGC measurement (stat.)"),
        Patch(facecolor="0.45", edgecolor="none", alpha=0.20,
              label="Point-to-point syst."),
        Line2D([], [], color="black", lw=2.0, label="Kroll prediction"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, handlelength=2.5)

    text = (
        "Kroll curves: virtual-photon frame\n"
        "RGC UL/LL points: lab-frame definition\n"
        "No frame transformation applied"
    )
    if x_index is not None:
        xlo, xhi = XB_BINS[x_index]
        text = rf"RGC: ${xlo:.2f} < x_B < {xhi:.2f}$" + "\n" + text
    # endif
    ax.text(0.02, 0.58, text, transform=ax.transAxes, va="top", ha="left",
            fontsize=10.5, linespacing=1.45)


def plot_one_xb_canvas(
    measurements: pd.DataFrame,
    predictions: pd.DataFrame,
    x_index: int,
    output_dir: Path,
) -> Path:
    subset = measurements.loc[measurements["x_index_plot"] == x_index].copy()
    kroll_xb = (0.2, 0.3, 0.4, 0.5)[x_index]
    pred = predictions.loc[np.isclose(predictions["xB"], kroll_xb)].copy()
    color = XB_COLORS[x_index]

    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.4), constrained_layout=False)
    axes_flat = axes.ravel()

    for i, observable in enumerate(OBSERVABLES):
        ax = axes_flat[i]
        draw_kroll_curve(ax, pred, observable, color="black")
        draw_measurements(ax, subset, observable, color=color, marker=XB_MARKERS[x_index])
        format_axis(ax, observable, show_xlabel=(i >= 3))
        lo, hi = robust_y_limits([subset], [pred], observable)
        ax.set_ylim(lo, hi)
        ax.text(0.04, 0.94, f"({chr(97 + i)})", transform=ax.transAxes,
                va="top", ha="left", fontweight="bold")
    # endfor

    info_panel(axes_flat[5], x_index=x_index)
    xlo, xhi = XB_BINS[x_index]
    mean_xb = subset["mean_xB"].mean()
    fig.suptitle(
        rf"RGC $e n\pi^+$ structure-function ratios: "
        rf"${xlo:.2f} < x_B < {xhi:.2f}$ ($\langle x_B\rangle={mean_xb:.3f}$)",
        y=0.985,
        fontsize=15,
    )
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.09, top=0.91,
                        wspace=0.34, hspace=0.28)

    path = output_dir / f"kroll_comparison_xB_bin_{x_index + 1}.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def plot_all_xb_canvas(
    measurements: pd.DataFrame,
    predictions: pd.DataFrame,
    output_dir: Path,
) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(12.8, 7.6), constrained_layout=False)
    axes_flat = axes.ravel()

    for i, observable in enumerate(OBSERVABLES):
        ax = axes_flat[i]
        measurement_subsets: list[pd.DataFrame] = []
        prediction_subsets: list[pd.DataFrame] = []

        for x_index, (color, marker) in enumerate(zip(XB_COLORS, XB_MARKERS)):
            subset = measurements.loc[measurements["x_index_plot"] == x_index].copy()
            kroll_xb = (0.2, 0.3, 0.4, 0.5)[x_index]
            pred = predictions.loc[np.isclose(predictions["xB"], kroll_xb)].copy()
            measurement_subsets.append(subset)
            prediction_subsets.append(pred)
            draw_kroll_curve(ax, pred, observable, color=color)
            draw_measurements(ax, subset, observable, color=color, marker=marker)
        # endfor

        format_axis(ax, observable, show_xlabel=(i >= 3))
        lo, hi = robust_y_limits(measurement_subsets, prediction_subsets, observable)
        ax.set_ylim(lo, hi)
        ax.text(0.04, 0.94, f"({chr(97 + i)})", transform=ax.transAxes,
                va="top", ha="left", fontweight="bold")
    # endfor

    legend_ax = axes_flat[5]
    legend_ax.axis("off")
    xb_handles = []
    for x_index, ((xlo, xhi), color, marker) in enumerate(zip(XB_BINS, XB_COLORS, XB_MARKERS)):
        xb_handles.append(
            Line2D([], [], color=color, lw=2.0, marker=marker, mfc="white",
                   mec=color, mew=1.2, ms=6,
                   label=rf"${xlo:.2f}<x_B<{xhi:.2f}$; Kroll $x_B={(0.2,0.3,0.4,0.5)[x_index]:.1f}$")
        )
    # endfor
    style_handles = [
        Line2D([], [], marker="o", linestyle="none", mfc="white", mec="black",
               mew=1.3, ms=6.5, label="RGC point: statistical error bar"),
        Patch(facecolor="0.45", edgecolor="none", alpha=0.20,
              label="RGC point-to-point systematic"),
        Line2D([], [], color="black", lw=2.0, label="Solid line: Kroll prediction"),
    ]
    legend1 = legend_ax.legend(handles=xb_handles, loc="upper left", frameon=False,
                               title=r"$x_B$ groups", title_fontsize=11)
    legend_ax.add_artist(legend1)
    legend_ax.legend(handles=style_handles, loc="lower left", frameon=False)
    legend_ax.text(
        0.02, 0.46,
        "Kroll: virtual-photon frame\nRGC UL/LL: lab-frame definition\nNo frame transformation applied",
        transform=legend_ax.transAxes, va="top", fontsize=10.3, linespacing=1.4,
    )

    fig.suptitle(r"RGC $e n\pi^+$ structure-function ratios and Kroll predictions",
                 y=0.985, fontsize=15)
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.09, top=0.91,
                        wspace=0.34, hspace=0.28)

    path = output_dir / "kroll_comparison_all_xB.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Final extraction CSV (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"PNG output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_style()

    input_csv = args.input_csv.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not input_csv.is_file():
        raise FileNotFoundError(
            f"Final extraction CSV not found: {input_csv}\n"
            "Run extract_structure_function_ratios.py first, or pass --input-csv."
        )
    # endif
    output_dir.mkdir(parents=True, exist_ok=True)

    measurements = validate_measurements(pd.read_csv(input_csv))
    predictions = kroll_frame()

    outputs: list[Path] = []
    for x_index in range(4):
        outputs.append(plot_one_xb_canvas(measurements, predictions, x_index, output_dir))
    # endfor
    outputs.append(plot_all_xb_canvas(measurements, predictions, output_dir))

    print("Kroll comparison plots complete.")
    print(f"  Input:  {input_csv}")
    print(f"  Output: {output_dir}")
    for path in outputs:
        print(f"    {path.name}")
    # endfor


if __name__ == "__main__":
    main()
