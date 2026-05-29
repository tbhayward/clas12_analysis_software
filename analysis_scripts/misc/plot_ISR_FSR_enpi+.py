#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ISR/FSR comparison and radiative-correction summary from fit-result text files.

This script reads two Mathematica-style fit-result text files:

    Baseline file:
        /u/home/thayward/tims_fits/RGC.txt

    ISR/FSR file:
        /u/home/thayward/tims_fits/RGC_ISR.txt

Each input line is expected to look like:

    enpiLowxBGEchi2FitsALUsinphi = {{x, y, sigma}, {x, y, sigma}, ...};

The script parses those lines into dictionaries of the form:

    data[bin_tag][series_name] = [[x, y, sigma], ...]

where, for example,

    bin_tag     = enpiLowxBGE
    series_name = ALUsinphi

It then generates, for each xB bin:
    1) A 2x3 comparison plot:
        subplot 1: F_LU^{sinphi}/F_UU
        subplot 2: F_UL^{sinphi}/F_UU
        subplot 3: F_UL^{sin2phi}/F_UU
        subplot 4: empty
        subplot 5: F_LL/F_UU
        subplot 6: F_LL^{cosphi}/F_UU

    2) A 2x3 signed-delta plot:
        Delta_rad = Baseline minus ISR&FSR

        Since the Baseline and ISR/FSR fits are made from the same underlying
        data, the delta points are plotted without statistical error bars.
        The two fit results are correlated, so the uncorrelated expression

            sigma_Delta = sqrt(sigma_Baseline^2 + sigma_ISRFSR^2)

        is intentionally not used.

    3) Per-bin CSV files with:
        -t' (GeV^2), Baseline value, Baseline sigma,
        ISR&FSR value, ISR&FSR sigma,
        Delta_rad_signed, abs_Delta_rad, sigma_rad

       where:
        Delta_rad_signed = Baseline - ISR&FSR
        abs_Delta_rad    = abs(Delta_rad_signed)
        sigma_rad        = sqrt(abs_Delta_rad^2 + EXTRA_RAD_SYS_ABS^2)

    4) A text radiative-correction summary in the block format expected by
       the downstream plotting script:
           -t'    Delta    sigma_rad

       The downstream script applies:
           A_Born = A_baseline - Delta

       and uses:
           sigma_rad = sqrt(abs(Delta)^2 + EXTRA_RAD_SYS_ABS^2)

    5) A flat machine-readable CSV:
           ISR_FSR_radiative_corrections.csv

    6) A sanity report checking:
        - missing series
        - length mismatches
        - identical series
        - zero deltas
        - x-grid mismatches

Important radiative-systematic convention:
    The signed radiative correction is still:

        Delta_rad = Baseline - ISR&FSR

    The assigned radiative systematic is now:

        sigma_rad = sqrt(abs(Delta_rad)^2 + 0.02^2)

    where 0.02 is an absolute 2 percent systematic uncertainty on the
    dimensionless asymmetry / structure-function ratio.

Dilution-factor panels and hard-coded dilution-factor arrays have been removed.
"""

import os
import re
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt

# =========================
# Default input/output paths
# =========================
DEFAULT_BASELINE_PATH = "/u/home/thayward/tims_fits/RGC.txt"
DEFAULT_ISRFSR_PATH = "/u/home/thayward/tims_fits/RGC_ISR.txt"
DEFAULT_OUT_DIR = os.path.join("output", "enpi+")

# =========================
# Radiative systematic convention
# =========================
# Reviewer-coordinated extra absolute radiative systematic uncertainty.
# This is interpreted as an absolute 2 percent uncertainty on the
# dimensionless fitted structure-function ratio.
#
# Final radiative systematic written to output:
#   sigma_rad = sqrt(abs(Delta_rad)^2 + EXTRA_RAD_SYS_ABS^2)
#
# The signed Delta_rad correction itself is NOT changed.
EXTRA_RAD_SYS_ABS = 0.02

# =========================
# Plot styling / constants
# =========================
COLORS = {
    "Baseline": "tab:blue",
    "ISR&FSR (data-driven)": "tab:red",
}

MARKERS = {
    "Baseline": "o",
    "ISR&FSR (data-driven)": "s",
}

CAPSIZE = 3
MS = 5.0
XLIM_T = (0.0, 1.30)
X_LABEL = r"$-t'\ (\mathrm{GeV}^{2})$"

YLIM_LU = (-0.40, 0.40)
YLIM_UL = (-0.40, 0.40)
YLIM_LL = (-1.00, 1.00)
YLIM_DELTA = (-0.30, 0.30)

TOP_PAD = 0.95
ATOL_ZERO = 1e-14

XB_LABELS = {
    "enpiGE":          r"$0.10 < x_{B} < 0.60$",
    "enpiLowxBGE":     r"$0.10 < x_{B} < 0.25$",
    "enpiMidLowxBGE":  r"$0.25 < x_{B} < 0.35$",
    "enpiMidHighxBGE": r"$0.35 < x_{B} < 0.45$",
    "enpiHighxBGE":    r"$0.45 < x_{B} < 0.60$",
}

SAVE_TAG = {
    "enpiGE": "enpi",
    "enpiLowxBGE": "enpiLowxB",
    "enpiMidLowxBGE": "enpiMidLowxB",
    "enpiMidHighxBGE": "enpiMidHighxB",
    "enpiHighxBGE": "enpiHighxB",
}

BIN_ORDER = [
    "enpiGE",
    "enpiLowxBGE",
    "enpiMidLowxBGE",
    "enpiMidHighxBGE",
    "enpiHighxBGE",
]

SERIES_ALL = [
    "ALUoffset",
    "AULoffset",
    "ALUsinphi",
    "AULsinphi",
    "AULsin2phi",
    "A_T_UL",
    "ALL",
    "ALLcosphi",
    "AUUcosphi",
    "AUUcos2phi",
    "A_T_LL",
]

# The five plotted physics series.
#
# Canvas layout:
#   subplot 1: BSA
#   subplot 2: TSA sin(phi)
#   subplot 3: TSA sin(2phi)
#   subplot 4: empty
#   subplot 5: DSA
#   subplot 6: DSA cos(phi)
SERIES_TO_PLOT = [
    ("ALUsinphi",  "ALUsin",  r"$F_{LU}^{\sin\phi}/F_{UU}$",   YLIM_LU),
    ("AULsinphi",  "AULsin",  r"$F_{UL}^{\sin\phi}/F_{UU}$",   YLIM_UL),
    ("AULsin2phi", "AULsin2", r"$F_{UL}^{\sin2\phi}/F_{UU}$",  YLIM_UL),
    ("ALL",        "ALLn0",   r"$F_{LL}/F_{UU}$",              YLIM_LL),
    ("ALLcosphi",  "ALLcos",  r"$F_{LL}^{\cos\phi}/F_{UU}$",   YLIM_LL),
]

# =========================
# Text parsing
# =========================
def read_text(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file does not exist: {path}")
    #endif

    with open(path, "r", encoding="utf-8") as f:
        return f.read()
    #endwith
#enddef


def parse_fit_text(path):
    """
    Parse Mathematica-style fit assignments from a text file.

    Expected assignment format:

        enpiLowxBGEchi2FitsALUsinphi = {{-1.1, 0.1, 0.01}, {...}};

    Returns:
        data[bin_tag][series_name] = [[x, y, sigma], ...]
    """
    text = read_text(path)

    data = {}

    assignment_pattern = re.compile(
        r"([A-Za-z0-9_]+chi2Fits[A-Za-z0-9_]+)\s*=\s*(\{\{.*?\}\})\s*;",
        re.DOTALL,
    )

    triple_pattern = re.compile(
        r"\{\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*,\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*,\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*"
        r"\}"
    )

    n_assignments = 0

    for match in assignment_pattern.finditer(text):
        full_name = match.group(1)
        array_text = match.group(2)

        split_token = "chi2Fits"

        if split_token not in full_name:
            continue
        #endif

        bin_tag, series_name = full_name.split(split_token, 1)

        triples = []

        for tmatch in triple_pattern.finditer(array_text):
            x = float(tmatch.group(1))
            y = float(tmatch.group(2))
            e = float(tmatch.group(3))
            triples.append([x, y, e])
        #endfor

        if len(triples) == 0:
            print(f"[WARN] No triples parsed for {full_name} in {path}.")
            continue
        #endif

        if bin_tag not in data:
            data[bin_tag] = {}
        #endif

        if series_name in data[bin_tag]:
            print(f"[WARN] Duplicate assignment for {bin_tag}:{series_name} in {path}; overwriting previous value.")
        #endif

        data[bin_tag][series_name] = triples
        n_assignments += 1
    #endfor

    if n_assignments == 0:
        raise RuntimeError(f"No fit assignments were parsed from file: {path}")
    #endif

    return data
#enddef


def summarize_loaded_data(label, data):
    n_bins = len(data)
    n_series = sum(len(v) for v in data.values())

    print(f"[OK] Loaded {label}: {n_bins} bins, {n_series} series")

    for b in sorted(data.keys()):
        available = sorted(data[b].keys())
        print(f"     {b}: {len(available)} series")
    #endfor
#enddef


# =========
# Helpers
# =========
def to_series(triples, negate_x=True, sort=True):
    """Convert list-of-[x,y,e] to aligned numpy arrays; optionally negate x to get -t'."""
    arr = np.array(triples, dtype=float)

    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"Expected an Nx3 array, got shape {arr.shape}")
    #endif

    x = -arr[:, 0] if negate_x else arr[:, 0]
    y = arr[:, 1]
    e = arr[:, 2]

    if sort:
        idx = np.argsort(x)
        x, y, e = x[idx], y[idx], e[idx]
    #endif

    return x, y, e
#enddef


def delta_rad(yb, yi):
    """Signed radiative correction shift: Baseline minus ISR&FSR."""
    d = yb - yi
    return d
#enddef


def compute_sigma_rad(delta_values):
    """
    Compute the radiative systematic uncertainty.

    The signed correction is:
        Delta_rad = Baseline - ISR&FSR

    The assigned uncertainty is:
        sigma_rad = sqrt(abs(Delta_rad)^2 + EXTRA_RAD_SYS_ABS^2)

    EXTRA_RAD_SYS_ABS is the reviewer-coordinated additional absolute
    2 percent systematic uncertainty used to cover unresolved radiative
    effects such as external radiation.
    """
    d = np.asarray(delta_values, dtype=float)
    return np.sqrt(np.abs(d) * np.abs(d) + EXTRA_RAD_SYS_ABS * EXTRA_RAD_SYS_ABS)
#enddef


def assert_same_length(bin_tag, series_name, xb, xi):
    if len(xb) != len(xi):
        raise ValueError(
            f"[ERROR] Bin={bin_tag} Series={series_name}: "
            f"mismatched counts (Baseline {len(xb)} vs ISR&FSR {len(xi)})."
        )
    #endif
#enddef


def get_bin_title(bin_tag):
    return r"$ep \rightarrow en\pi^{+}$ -- " + XB_LABELS.get(bin_tag, bin_tag)
#enddef


def get_existing_bin_tags(baseline, isrfsr):
    tags = []

    for b in BIN_ORDER:
        if b in baseline and b in isrfsr:
            tags.append(b)
        #endif
    #endfor

    extra = sorted((set(baseline.keys()) & set(isrfsr.keys())) - set(tags))

    for b in extra:
        tags.append(b)
    #endfor

    return tags
#enddef


# =========
# Plotters
# =========
def legend_on(ax, loc="upper right"):
    h, l = ax.get_legend_handles_labels()

    if h:
        leg = ax.legend(h, l, frameon=True, edgecolor="black", fontsize=11, loc=loc)
        leg.get_frame().set_alpha(0.9)
    #endif
#enddef


def draw_compare(ax, baseline, isrfsr, bin_tag, series_name, ylabel, ylim=None):
    sources = [
        ("Baseline", baseline),
        ("ISR&FSR (data-driven)", isrfsr),
    ]

    for lab, src in sources:
        triples = src.get(bin_tag, {}).get(series_name)

        if not triples:
            continue
        #endif

        x, y, e = to_series(triples)

        ax.errorbar(
            x,
            y,
            yerr=e,
            fmt=MARKERS[lab],
            color=COLORS[lab],
            ecolor=COLORS[lab],
            capsize=CAPSIZE,
            markersize=MS,
            linestyle="None",
            label=lab,
        )
    #endfor

    if ylim is not None:
        ax.set(ylim=ylim)
    #endif

    ax.set(xlim=XLIM_T, xlabel=X_LABEL, ylabel=ylabel)
    ax.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ax.grid(True, linestyle="--", alpha=0.6)
#enddef


def draw_delta_indexed(ax, baseline, isrfsr, bin_tag, series_name, ylabel):
    b = baseline.get(bin_tag, {}).get(series_name)
    i = isrfsr.get(bin_tag, {}).get(series_name)

    if not b or not i:
        return
    #endif

    xb, yb, eb = to_series(b)
    xi, yi, ei = to_series(i)

    assert_same_length(bin_tag, series_name, xb, xi)

    if not np.allclose(xb, xi, rtol=0, atol=1e-6):
        print(f"[WARN] x mismatch in {bin_tag}:{series_name} (index-paired).")
    #endif

    yd = delta_rad(yb, yi)

    ax.plot(
        xb,
        yd,
        "o",
        color="k",
        markersize=MS,
        linestyle="None",
    )

    ax.set(xlim=XLIM_T, ylim=YLIM_DELTA, xlabel=X_LABEL, ylabel=ylabel)
    ax.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ax.grid(True, linestyle="--", alpha=0.6)
#enddef


def plot_bin(baseline, isrfsr, bin_tag, out_dir):
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.6))
    ax_lu, ax_ul1, ax_ul2, ax_empty, ax_ll0, ax_ll1 = axes.flat

    draw_compare(
        ax_lu,
        baseline,
        isrfsr,
        bin_tag,
        "ALUsinphi",
        r"$F_{LU}^{\sin\phi}/F_{UU}$",
        YLIM_LU,
    )

    draw_compare(
        ax_ul1,
        baseline,
        isrfsr,
        bin_tag,
        "AULsinphi",
        r"$F_{UL}^{\sin\phi}/F_{UU}$",
        YLIM_UL,
    )

    draw_compare(
        ax_ul2,
        baseline,
        isrfsr,
        bin_tag,
        "AULsin2phi",
        r"$F_{UL}^{\sin2\phi}/F_{UU}$",
        YLIM_UL,
    )

    ax_empty.axis("off")

    draw_compare(
        ax_ll0,
        baseline,
        isrfsr,
        bin_tag,
        "ALL",
        r"$F_{LL}/F_{UU}$",
        YLIM_LL,
    )

    draw_compare(
        ax_ll1,
        baseline,
        isrfsr,
        bin_tag,
        "ALLcosphi",
        r"$F_{LL}^{\cos\phi}/F_{UU}$",
        YLIM_LL,
    )

    legend_on(ax_ul1, loc="upper right")
    legend_on(ax_ll1, loc="upper right")

    plt.suptitle(get_bin_title(bin_tag), fontsize=16, y=0.97)
    plt.tight_layout(rect=[0, 0, 1, TOP_PAD])

    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"ISR_FSR_comparison_{SAVE_TAG.get(bin_tag, bin_tag)}.pdf")

    plt.savefig(path)
    plt.close(fig)

    print(f"[OK] Saved: {path}")
#enddef


def plot_bin_delta(baseline, isrfsr, bin_tag, out_dir):
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.6))
    ax_lu, ax_ul1, ax_ul2, ax_empty, ax_ll0, ax_ll1 = axes.flat

    draw_delta_indexed(
        ax_lu,
        baseline,
        isrfsr,
        bin_tag,
        "ALUsinphi",
        r"$\Delta_{\mathrm{rad}}(F_{LU}^{\sin\phi}/F_{UU})$",
    )

    draw_delta_indexed(
        ax_ul1,
        baseline,
        isrfsr,
        bin_tag,
        "AULsinphi",
        r"$\Delta_{\mathrm{rad}}(F_{UL}^{\sin\phi}/F_{UU})$",
    )

    draw_delta_indexed(
        ax_ul2,
        baseline,
        isrfsr,
        bin_tag,
        "AULsin2phi",
        r"$\Delta_{\mathrm{rad}}(F_{UL}^{\sin2\phi}/F_{UU})$",
    )

    ax_empty.axis("off")

    draw_delta_indexed(
        ax_ll0,
        baseline,
        isrfsr,
        bin_tag,
        "ALL",
        r"$\Delta_{\mathrm{rad}}(F_{LL}/F_{UU})$",
    )

    draw_delta_indexed(
        ax_ll1,
        baseline,
        isrfsr,
        bin_tag,
        "ALLcosphi",
        r"$\Delta_{\mathrm{rad}}(F_{LL}^{\cos\phi}/F_{UU})$",
    )

    plt.suptitle(get_bin_title(bin_tag), fontsize=16, y=0.97)
    plt.tight_layout(rect=[0, 0, 1, TOP_PAD])

    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"ISR_FSR_delta_{SAVE_TAG.get(bin_tag, bin_tag)}.pdf")

    plt.savefig(path)
    plt.close(fig)

    print(f"[OK] Saved: {path}")
#enddef


# =================
# Text/CSV outputs
# =================
def collect_delta_sf(baseline, isrfsr, bin_tag, series_name):
    b = baseline.get(bin_tag, {}).get(series_name)
    i = isrfsr.get(bin_tag, {}).get(series_name)

    if not b or not i:
        return None
    #endif

    xb, yb, eb = to_series(b)
    xi, yi, ei = to_series(i)

    assert_same_length(bin_tag, series_name, xb, xi)

    if not np.allclose(xb, xi, rtol=0, atol=1e-6):
        print(f"[WARN] x mismatch in {bin_tag}:{series_name} while collecting delta.")
    #endif

    d = delta_rad(yb, yi)
    abs_d = np.abs(d)
    s_rad = compute_sigma_rad(d)

    return xb, yb, eb, yi, ei, d, abs_d, s_rad
#enddef


def write_delta_summary(out_dir, baseline, isrfsr, bin_tags):
    """
    Write the block-format summary expected by the downstream plotting script.

    Format:
        Bin: <bin_tag> ...
        Series: <LaTeX label>
             -t'       Delta    sigma_rad
          0.12345   signed_D   sqrt(abs(signed_D)^2 + 0.02^2)

    The downstream script applies:
        y <- y - Delta

    and uses:
        sigma_rad = sqrt(abs(Delta)^2 + EXTRA_RAD_SYS_ABS^2)
    """
    path = os.path.join(out_dir, "ISR_FSR_delta_summary.txt")

    with open(path, "w", encoding="utf-8") as f:
        f.write("Radiative correction summary: Delta_rad = Baseline minus ISR&FSR\n")
        f.write("Units: -t' in GeV^2\n")
        f.write("\n")
        f.write("Correction convention:\n")
        f.write("  The downstream correction should use A_Born = A_baseline - Delta_rad.\n")
        f.write("  Delta_rad is signed and is listed in the column named Delta.\n")
        f.write("\n")
        f.write("Uncertainty convention:\n")
        f.write("  Baseline and ISR&FSR are fits to the same underlying data.\n")
        f.write("  Therefore the uncorrelated expression sqrt(sigma_Baseline^2 + sigma_ISRFSR^2) is not used.\n")
        f.write("  The assigned radiative systematic is sigma_rad = sqrt(abs(Delta_rad)^2 + extra_rad_sys^2).\n")
        f.write("  extra_rad_sys = {:.6f}, an absolute 2 percent bin-by-bin systematic uncertainty.\n".format(EXTRA_RAD_SYS_ABS))
        f.write("  This extra term covers unresolved/undetermined radiative effects, including external radiation.\n")
        f.write("=" * 88 + "\n")

        col_t = "-t'"
        col_d = "Delta"
        col_s = "sigma_rad"

        for b in bin_tags:
            f.write(f"\nBin: {b}    x_B range: {XB_LABELS.get(b, '')}\n")
            f.write("-" * 88 + "\n")

            for series_name, series_key, label, _ in SERIES_TO_PLOT:
                res = collect_delta_sf(baseline, isrfsr, b, series_name)

                if res is None:
                    f.write(f"Series: {label}\n")
                    f.write("MISSING in Baseline and/or ISR&FSR input.\n\n")
                    continue
                #endif

                x, _, _, _, _, d, abs_d, s_rad = res

                f.write(f"Series: {label}\n")
                f.write("{:>8}    {:>12}    {:>12}\n".format(col_t, col_d, col_s))

                for xi, di_, si_ in zip(x, d, s_rad):
                    f.write(f"{xi:8.5f}    {di_:12.6f}    {si_:12.6f}\n")
                #endfor

                f.write("\n")
            #endfor
        #endfor
    #endwith

    print(f"[OK] Wrote radiative correction summary: {path}")
#enddef


def write_machine_readable_corrections_csv(out_dir, baseline, isrfsr, bin_tags):
    """
    Write a flat CSV that is easier to feed into future scripts.

    The essential columns are:
        delta_rad_signed = Baseline - ISR&FSR
        abs_delta_rad    = abs(delta_rad_signed)
        sigma_rad        = sqrt(abs_delta_rad^2 + EXTRA_RAD_SYS_ABS^2)

    This file is intentionally redundant with ISR_FSR_delta_summary.txt, but is
    more convenient for future programmatic use.
    """
    path = os.path.join(out_dir, "ISR_FSR_radiative_corrections.csv")

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)

        w.writerow([
            "bin_tag",
            "xB_range",
            "series_key",
            "series_label",
            "-t' (GeV^2)",
            "baseline",
            "sigma_baseline",
            "isrfsr",
            "sigma_isrfsr",
            "delta_rad_signed",
            "abs_delta_rad",
            "extra_rad_sys_abs",
            "sigma_rad",
        ])

        for b in bin_tags:
            for series_name, series_key, label, _ in SERIES_TO_PLOT:
                res = collect_delta_sf(baseline, isrfsr, b, series_name)

                if res is None:
                    continue
                #endif

                x, yb, eb, yi, ei, d, abs_d, s_rad = res

                for row in zip(x, yb, eb, yi, ei, d, abs_d, s_rad):
                    w.writerow([
                        b,
                        XB_LABELS.get(b, ""),
                        series_key,
                        label,
                        f"{row[0]:.6f}",
                        f"{row[1]:.9g}",
                        f"{row[2]:.9g}",
                        f"{row[3]:.9g}",
                        f"{row[4]:.9g}",
                        f"{row[5]:.9g}",
                        f"{row[6]:.9g}",
                        f"{EXTRA_RAD_SYS_ABS:.9g}",
                        f"{row[7]:.9g}",
                    ])
                #endfor
            #endfor
        #endfor
    #endwith

    print(f"[OK] Wrote machine-readable correction CSV: {path}")
#enddef


def write_bin_csvs(out_dir, baseline, isrfsr, bin_tag):
    """Per-series CSV with x, baseline, isrfsr, signed delta, and sigma_rad."""
    bin_dir = os.path.join(out_dir, f"csv_{SAVE_TAG.get(bin_tag, bin_tag)}")
    os.makedirs(bin_dir, exist_ok=True)

    for series_name, series_key, label, _ in SERIES_TO_PLOT:
        res = collect_delta_sf(baseline, isrfsr, bin_tag, series_name)

        if res is None:
            print(f"[WARN] CSV skipped for {bin_tag}:{series_name} because it is missing in one input.")
            continue
        #endif

        x, yb, eb, yi, ei, d, abs_d, s_rad = res
        fname = series_name + ".csv"

        with open(os.path.join(bin_dir, fname), "w", newline="") as fh:
            w = csv.writer(fh)

            w.writerow([
                "-t' (GeV^2)",
                f"Baseline {label}",
                "sigma(B)",
                f"ISR&FSR {label}",
                "sigma(I)",
                "Delta_rad_signed",
                "abs_Delta_rad",
                "extra_rad_sys_abs",
                "sigma_rad",
            ])

            for row in zip(x, yb, eb, yi, ei, d, abs_d, s_rad):
                w.writerow([
                    f"{row[0]:.6f}",
                    f"{row[1]:.6g}",
                    f"{row[2]:.6g}",
                    f"{row[3]:.6g}",
                    f"{row[4]:.6g}",
                    f"{row[5]:.6g}",
                    f"{row[6]:.6g}",
                    f"{EXTRA_RAD_SYS_ABS:.6g}",
                    f"{row[7]:.6g}",
                ])
            #endfor
        #endwith
    #endfor
#enddef


def write_abs_shift_summary(out_dir, baseline, isrfsr, bin_tags):
    """
    Write a compact text summary of absolute shift magnitudes and final
    radiative systematic magnitudes:
        [min, mean, max]
    """
    path = os.path.join(out_dir, "ISR_FSR_abs_shift_summary.txt")

    with open(path, "w", encoding="utf-8") as f:
        f.write("Absolute radiative shift summary\n")
        f.write("For abs(Delta_rad), each entry is [min, mean, max].\n")
        f.write("For sigma_rad, each entry is [min, mean, max].\n")
        f.write("sigma_rad = sqrt(abs(Delta_rad)^2 + extra_rad_sys^2).\n")
        f.write("extra_rad_sys = {:.6f}, an absolute 2 percent bin-by-bin systematic uncertainty.\n".format(EXTRA_RAD_SYS_ABS))
        f.write("=" * 88 + "\n")

        for b in bin_tags:
            f.write(f"\nBin: {b}    x_B range: {XB_LABELS.get(b, '')}\n")
            f.write("-" * 88 + "\n")

            for series_name, series_key, label, _ in SERIES_TO_PLOT:
                res = collect_delta_sf(baseline, isrfsr, b, series_name)

                if res is None:
                    continue
                #endif

                _, _, _, _, _, _, abs_d, s_rad = res
                f.write(
                    f"{series_key:8s}  {label:40s}  "
                    f"absDelta=[{np.min(abs_d):.6f}, {np.mean(abs_d):.6f}, {np.max(abs_d):.6f}]  "
                    f"sigmaRad=[{np.min(s_rad):.6f}, {np.mean(s_rad):.6f}, {np.max(s_rad):.6f}]\n"
                )
            #endfor
        #endfor
    #endwith

    print(f"[OK] Wrote absolute shift summary: {path}")
#enddef


# ==========================
# Sanity / debug report
# ==========================
def sanity_lines_for_bin(baseline, isrfsr, bin_tag):
    lines = [f"== Bin: {bin_tag} ({XB_LABELS.get(bin_tag, '')})"]

    check_series = [
        "ALUsinphi",
        "AULsinphi",
        "AULsin2phi",
        "ALL",
        "ALLcosphi",
    ]

    for s in check_series:
        b = baseline.get(bin_tag, {}).get(s)
        i = isrfsr.get(bin_tag, {}).get(s)

        if not b or not i:
            missing = []

            if not b:
                missing.append("BASELINE")
            #endif

            if not i:
                missing.append("ISRFSR")
            #endif

            lines.append(f"  {s}: MISSING in {', '.join(missing)}")
            continue
        #endif

        xb, yb, eb = to_series(b)
        xi, yi, ei = to_series(i)

        same_len = len(yb) == len(yi)
        exact_equal = same_len and np.array_equal(yb, yi)
        close_equal = same_len and np.allclose(yb, yi)

        if same_len:
            d = delta_rad(yb, yi)
            zero_delta = d.size > 0 and np.allclose(d, 0.0, atol=ATOL_ZERO)
            max_abs_d = float(np.max(np.abs(d))) if d.size else float("nan")
            mean_abs_d = float(np.mean(np.abs(d))) if d.size else float("nan")
            s_rad = compute_sigma_rad(d)
            max_s_rad = float(np.max(s_rad)) if s_rad.size else float("nan")
            mean_s_rad = float(np.mean(s_rad)) if s_rad.size else float("nan")
        else:
            zero_delta = False
            max_abs_d = float("nan")
            mean_abs_d = float("nan")
            max_s_rad = float("nan")
            mean_s_rad = float("nan")
        #endif

        lines.append(
            f"  {s}: lenB={len(yb)} lenI={len(yi)} "
            f"exact_equal={exact_equal} close_equal={close_equal} "
            f"Delta_all_zero={zero_delta} "
            f"mean|Delta|={mean_abs_d:.6g} max|Delta|={max_abs_d:.6g} "
            f"meanSigmaRad={mean_s_rad:.6g} maxSigmaRad={max_s_rad:.6g}"
        )

        if same_len and not np.allclose(xb, xi, atol=1e-6):
            lines.append("     NOTE: x grids differ (still paired by index).")
        #endif

        if exact_equal:
            lines.append("     >>> WARNING: y_Baseline and y_ISRFSR are IDENTICAL.")
        #endif
    #endfor

    return lines
#enddef


def write_sanity_report(out_dir, baseline, isrfsr, bin_tags):
    path = os.path.join(out_dir, "ISR_FSR_sanity_report.txt")

    with open(path, "w", encoding="utf-8") as f:
        f.write("Sanity report (index-based Delta_rad pairing)\n")
        f.write("Flags identical series / zero-Delta after pairing / x-grid mismatches.\n")
        f.write("Delta_rad is computed as Baseline minus ISR&FSR central values.\n")
        f.write("sigma_rad is assigned as sqrt(abs(Delta_rad)^2 + extra_rad_sys^2).\n")
        f.write("extra_rad_sys = {:.6f}, an absolute 2 percent bin-by-bin systematic uncertainty.\n\n".format(EXTRA_RAD_SYS_ABS))

        for b in bin_tags:
            for line in sanity_lines_for_bin(baseline, isrfsr, b):
                f.write(line + "\n")
            #endfor

            f.write("\n")
        #endfor
    #endwith

    print(f"[OK] Wrote sanity report: {path}")
#enddef


# ==========================
# Argument parsing
# ==========================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Make ISR/FSR comparison and radiative-correction outputs from two fit-result text files."
    )

    parser.add_argument(
        "--baseline",
        default=DEFAULT_BASELINE_PATH,
        help=f"Baseline fit-result text file. Default: {DEFAULT_BASELINE_PATH}",
    )

    parser.add_argument(
        "--isrfsr",
        default=DEFAULT_ISRFSR_PATH,
        help=f"ISR/FSR fit-result text file. Default: {DEFAULT_ISRFSR_PATH}",
    )

    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUT_DIR}",
    )

    return parser.parse_args()
#enddef


# =====
# Main
# =====
def main():
    args = parse_args()

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    print("[INFO] Radiative systematic convention:")
    print("[INFO]   Delta_rad = Baseline - ISR&FSR")
    print("[INFO]   Downstream central-value correction: A_Born = A_baseline - Delta_rad")
    print("[INFO]   sigma_rad = sqrt(abs(Delta_rad)^2 + extra_rad_sys^2)")
    print("[INFO]   extra_rad_sys = {:.6f} absolute, i.e. 2 percent bin-by-bin systematic".format(EXTRA_RAD_SYS_ABS))
    print("[INFO]   The signed Delta_rad values are unchanged by this extra uncertainty.")
    print("[INFO]   Only the sigma_rad output column is enlarged by the quadrature term.")
    print("")

    baseline = parse_fit_text(args.baseline)
    isrfsr = parse_fit_text(args.isrfsr)

    summarize_loaded_data("Baseline", baseline)
    summarize_loaded_data("ISR&FSR", isrfsr)

    bin_tags = get_existing_bin_tags(baseline, isrfsr)

    if len(bin_tags) == 0:
        raise RuntimeError("No common bin tags found between the two input files.")
    #endif

    print("[OK] Common bins to process:")

    for b in bin_tags:
        print(f"     {b}")
    #endfor

    write_sanity_report(out_dir, baseline, isrfsr, bin_tags)

    for b in bin_tags:
        plot_bin(baseline, isrfsr, b, out_dir)
        plot_bin_delta(baseline, isrfsr, b, out_dir)
        write_bin_csvs(out_dir, baseline, isrfsr, b)
    #endfor

    write_delta_summary(out_dir, baseline, isrfsr, bin_tags)
    write_machine_readable_corrections_csv(out_dir, baseline, isrfsr, bin_tags)
    write_abs_shift_summary(out_dir, baseline, isrfsr, bin_tags)

    print("")
    print("[OK] Completed ISR/FSR radiative correction outputs.")
    print("[OK] Reminder: sigma_rad now includes the extra 0.02 absolute systematic in quadrature.")
#enddef


if __name__ == "__main__":
    main()
#endif