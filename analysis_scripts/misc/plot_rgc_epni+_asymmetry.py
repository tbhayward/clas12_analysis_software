#!/usr/bin/env python3
import sys
import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------
# Usage:
#   python plot_enpi_from_texts.py Su22.txt Fa22.txt Sp23.txt Combined.txt enpiHightGE
#
# Saves:
#   output/enpi+/rgc_<PREFIX>_AllPeriods.pdf
#   output/enpi+/rgc_<PREFIX>_CombinedOnly.pdf
# ---------------------------------------------------------------------

# Colors for each run period (keep existing palette)
COLORS = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green",
    "Combined": "black",
}

# Marker/legend aesthetics (consistent with your earlier figures)
MARKER = "o"
CAPSIZE = 3
LABEL_FONTSIZE = 13

# ----------- parsing -------------------------------------------------
_assign_re = re.compile(r'([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};', re.DOTALL)
_triple_re = re.compile(r'\{([^{}]+)\}')

def parse_asym_file(path):
    """
    Parse a text file of the form:
      NAME = {{mean, val, err}, {mean, val, err}, ...};
    Return dict: name -> np.array(N,3): columns (mean, value, err)
    """
    with open(path, "r") as f:
        text = f.read()

    out = {}
    for m in _assign_re.finditer(text):
        name = m.group(1)
        content = m.group(2)
        triples = []
        for t in _triple_re.findall(content):
            parts = [p.strip() for p in t.split(",")]
            if len(parts) != 3:
                continue
            try:
                a = float(parts[0])
                b = float(parts[1])
                c = float(parts[2])
                # keep NaNs if present; filter later
                triples.append((a, b, c))
            except ValueError:
                continue
        if triples:
            out[name] = np.array(triples, dtype=float)
    return out

def get_series(dct, key):
    """
    Return dict(x,y,yerr) from parsed dict if key exists & non-empty; else None.
    Filters rows with non-finite y or non-positive/NaN yerr.
    """
    if key not in dct:
        return None
    arr = np.array(dct[key], dtype=float)
    if arr.size == 0:
        return None
    # filter invalid points (err <= 0 or NaN; y NaN; x NaN)
    x, y, e = arr[:,0], arr[:,1], arr[:,2]
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    if not np.any(mask):
        return None
    x, y, e = x[mask], y[mask], e[mask]
    return {"x": x, "y": y, "yerr": e}

# ----------- y-limits helper ----------------------------------------
def auto_ylim(datasets, pad_frac=0.12, min_span=1e-3):
    """
    datasets: iterable of dicts like {"y":..., "yerr":...} (None allowed)
    Returns (ymin, ymax) computed from min(y - e), max(y + e) with padding.
    """
    ymin = +np.inf
    ymax = -np.inf
    any_data = False
    for d in datasets:
        if d is None: 
            continue
        y, e = d["y"], d["yerr"]
        ylo = np.min(y - e) if y.size else +np.inf
        yhi = np.max(y + e) if y.size else -np.inf
        if np.isfinite(ylo): ymin = min(ymin, ylo)
        if np.isfinite(yhi): ymax = max(ymax, yhi)
        any_data = True
    if not any_data or not np.isfinite(ymin) or not np.isfinite(ymax):
        return (-1.0, 1.0)  # fallback
    # ensure some span
    span = max(ymax - ymin, min_span)
    pad = pad_frac * span
    return (ymin - pad, ymax + pad)

# ----------- building period dicts ----------------------------------
def build_period_dict(parsed, prefix):
    """
    For a parsed file dict and a prefix like 'enpiHightGE',
    return a mapping for the series needed per channel.
    """
    key = lambda suffix: f"{prefix}chi2Fits{suffix}"
    return {
        # UU (unpolarized denominator)
        "UUcos":   get_series(parsed, key("AUUcosphi")),
        "UUcos2":  get_series(parsed, key("AUUcos2phi")),
        # BSA
        "ALUsin":  get_series(parsed, key("ALUsinphi")),
        # TSA
        "AULsin":  get_series(parsed, key("AULsinphi")),
        "AULsin2": get_series(parsed, key("AULsin2phi")),
        # DSA
        "ALLn0":   get_series(parsed, key("ALL")),          # n=0 term
        "ALLcos":  get_series(parsed, key("ALLcosphi")),    # n=1 (cos φ)
    }

# ----------- plotting helpers ---------------------------------------
def plot_all_periods(p_su22, p_fa22, p_sp23, prefix, out_dir):
    plt.figure(figsize=(12, 9))
    plt.suptitle(
        r"$ep \rightarrow en\pi^{+}$, $0.07<|t|<0.7$, $z>0.55$, $y<0.65$, $0.75<M_{x}^{2}<1.05\ \mathrm{GeV}^{2}$",
        fontsize=16, y=0.97
    )

    # Axes
    axUU  = plt.subplot(2,2,1)  # UU (cos, cos2)
    axLU  = plt.subplot(2,2,2)  # BSA
    axUL  = plt.subplot(2,2,3)  # TSA
    axLL  = plt.subplot(2,2,4)  # DSA

    # helper to draw one period with styles
    def draw_period(ax, series, color, label, harmonic_style=None):
        """
        series: dict with keys present or None.
        harmonic_style:
          - for 2-trace plots: ("open_key", "filled_key")
        """
        if harmonic_style is None:
            # single series expected under 'single'
            s = series
            if s is not None:
                ax.errorbar(s["x"], s["y"], yerr=s["yerr"],
                            fmt=MARKER, color=color, ecolor=color,
                            mfc=color, mec=color, capsize=CAPSIZE, label=label)
        else:
            open_key, filled_key = harmonic_style
            s_open  = series.get(open_key, None)
            s_fill  = series.get(filled_key, None)
            if s_open is not None:
                ax.errorbar(s_open["x"], s_open["y"], yerr=s_open["yerr"],
                            fmt=MARKER, mfc="none", mec=color, ecolor=color,
                            capsize=CAPSIZE, label=label)  # open marker
            if s_fill is not None:
                ax.errorbar(s_fill["x"], s_fill["y"], yerr=s_fill["yerr"],
                            fmt=MARKER, color=color, ecolor=color,
                            capsize=CAPSIZE)  # filled marker; no extra label

    # ---- Top-left: UU (cos nφ) - n=1 open, n=2 filled (like TSA style) ----
    for lab, pdata in (("Su22", p_su22), ("Fa22", p_fa22), ("Sp23", p_sp23)):
        draw_period(axUU, {"UUcos": pdata["UUcos"], "UUcos2": pdata["UUcos2"]},
                    COLORS[lab], lab, harmonic_style=("UUcos", "UUcos2"))
    axUU.set(xlim=(0, 0.7), xlabel=r"$x_{B}$", ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ylo,yhi = auto_ylim([p_su22["UUcos"], p_su22["UUcos2"],
                         p_fa22["UUcos"], p_fa22["UUcos2"],
                         p_sp23["UUcos"], p_sp23["UUcos2"]])
    axUU.set_ylim(ylo, yhi)
    axUU.grid(True, linestyle="--", alpha=0.6)
    # legends
    harmUU = axUU.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUU.add_artist(harmUU)
    legUU = axUU.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, color=COLORS["Su22"], linestyle='', label='Su22'),
            Line2D([0],[0], marker=MARKER, color=COLORS["Fa22"], linestyle='', label='Fa22'),
            Line2D([0],[0], marker=MARKER, color=COLORS["Sp23"], linestyle='', label='Sp23'),
        ],
        title="Run Period", frameon=True, edgecolor="black", loc="lower right",
        fontsize=11, title_fontsize=12
    )
    legUU.get_frame().set_alpha(0.9)

    # ---- Top-right: BSA (ALU sinφ) ----
    for lab, pdata in (("Su22", p_su22), ("Fa22", p_fa22), ("Sp23", p_sp23)):
        s = pdata["ALUsin"]
        if s is not None:
            axLU.errorbar(s["x"], s["y"], yerr=s["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=lab)
    axLU.set(xlim=(0, 0.7), xlabel=r"$x_{B}$", ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ylo,yhi = auto_ylim([p_su22["ALUsin"], p_fa22["ALUsin"], p_sp23["ALUsin"]])
    axLU.set_ylim(ylo, yhi)
    axLU.grid(True, linestyle="--", alpha=0.6)
    legLU = axLU.legend(title="Run Period", frameon=True, edgecolor="black",
                        fontsize=11, title_fontsize=12)
    legLU.get_frame().set_alpha(0.9)

    # ---- Bottom-left: TSA (AUL sinφ open, sin2φ filled) ----
    for lab, pdata in (("Su22", p_su22), ("Fa22", p_fa22), ("Sp23", p_sp23)):
        # n=1 open
        s1 = pdata["AULsin"]
        if s1 is not None:
            axUL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt=MARKER, mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1")
        # n=2 filled
        s2 = pdata["AULsin2"]
        if s2 is not None:
            axUL.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=2")
    axUL.set(xlim=(0, 0.7), xlabel=r"$x_{B}$", ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ylo,yhi = auto_ylim([p_su22["AULsin"], p_su22["AULsin2"],
                         p_fa22["AULsin"], p_fa22["AULsin2"],
                         p_sp23["AULsin"], p_sp23["AULsin2"]])
    axUL.set_ylim(ylo, yhi)
    axUL.grid(True, linestyle="--", alpha=0.6)
    # harmonic legend
    harmUL = axUL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2')
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(harmUL)
    # run-period legend
    legUL = axUL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, color=COLORS["Su22"], linestyle='', label='Su22'),
            Line2D([0],[0], marker=MARKER, color=COLORS["Fa22"], linestyle='', label='Fa22'),
            Line2D([0],[0], marker=MARKER, color=COLORS["Sp23"], linestyle='', label='Sp23')
        ],
        title="Run Period", frameon=True, edgecolor="black", loc="lower right",
        fontsize=11, title_fontsize=12
    )
    legUL.get_frame().set_alpha(0.9)

    # ---- Bottom-right: DSA (ALL n=0 open, cosφ filled) ----
    for lab, pdata in (("Su22", p_su22), ("Fa22", p_fa22), ("Sp23", p_sp23)):
        s0 = pdata["ALLn0"]
        s1 = pdata["ALLcos"]
        if s0 is not None:
            axLL.errorbar(s0["x"], s0["y"], yerr=s0["yerr"],
                          fmt=MARKER, mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=0")
        if s1 is not None:
            axLL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1")
    axLL.set(xlim=(0, 0.7), xlabel=r"$x_{B}$", ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ylo,yhi = auto_ylim([p_su22["ALLn0"], p_su22["ALLcos"],
                         p_fa22["ALLn0"], p_fa22["ALLcos"],
                         p_sp23["ALLn0"], p_sp23["ALLcos"]])
    axLL.set_ylim(ylo, yhi)
    axLL.grid(True, linestyle="--", alpha=0.6)
    harmLL = axLL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=1')
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axLL.add_artist(harmLL)
    legLL = axLL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, color=COLORS["Su22"], linestyle='', label='Su22'),
            Line2D([0],[0], marker=MARKER, color=COLORS["Fa22"], linestyle='', label='Fa22'),
            Line2D([0],[0], marker=MARKER, color=COLORS["Sp23"], linestyle='', label='Sp23')
        ],
        title="Run Period", frameon=True, edgecolor="black", loc="lower right",
        fontsize=11, title_fontsize=12
    )
    legLL.get_frame().set_alpha(0.9)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{prefix}_AllPeriods.pdf")
    plt.savefig(out_path)
    print(f"Saved all-periods figure: {out_path}")

def plot_combined_only(p_comb, prefix, out_dir):
    plt.figure(figsize=(12, 9))
    plt.suptitle(
        r"$ep \rightarrow en\pi^{+}$, $0.07<|t|<0.7$, $z>0.55$, $y<0.65$, $0.75<M_{x}^{2}<1.05\ \mathrm{GeV}^{2}$",
        fontsize=16, y=0.97
    )

    axUU  = plt.subplot(2,2,1)
    axLU  = plt.subplot(2,2,2)
    axUL  = plt.subplot(2,2,3)
    axLL  = plt.subplot(2,2,4)

    black = COLORS["Combined"]

    # UU: n=1 open, n=2 filled
    if p_comb["UUcos"] is not None:
        axUU.errorbar(p_comb["UUcos"]["x"], p_comb["UUcos"]["y"], yerr=p_comb["UUcos"]["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE)
    if p_comb["UUcos2"] is not None:
        axUU.errorbar(p_comb["UUcos2"]["x"], p_comb["UUcos2"]["y"], yerr=p_comb["UUcos2"]["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    axUU.set(xlim=(0, 0.7), xlabel=r"$x_{B}$", ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ylo,yhi = auto_ylim([p_comb["UUcos"], p_comb["UUcos2"]])
    axUU.set_ylim(ylo, yhi)
    axUU.grid(True, linestyle="--", alpha=0.6)
    harmUU = axUU.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUU.add_artist(harmUU)

    # BSA
    if p_comb["ALUsin"] is not None:
        s = p_comb["ALUsin"]
        axLU.errorbar(s["x"], s["y"], yerr=s["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    axLU.set(xlim=(0, 0.7), xlabel=r"$x_{B}$", ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ylo,yhi = auto_ylim([p_comb["ALUsin"]])
    axLU.set_ylim(ylo, yhi)
    axLU.grid(True, linestyle="--", alpha=0.6)

    # TSA
    if p_comb["AULsin"] is not None:
        s1 = p_comb["AULsin"]
        axUL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE)
    if p_comb["AULsin2"] is not None:
        s2 = p_comb["AULsin2"]
        axUL.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    axUL.set(xlim=(0, 0.7), xlabel=r"$x_{B}$", ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ylo,yhi = auto_ylim([p_comb["AULsin"], p_comb["AULsin2"]])
    axUL.set_ylim(ylo, yhi)
    axUL.grid(True, linestyle="--", alpha=0.6)
    harmUL = axUL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2'),
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(harmUL)

    # DSA
    if p_comb["ALLn0"] is not None:
        s0 = p_comb["ALLn0"]
        axLL.errorbar(s0["x"], s0["y"], yerr=s0["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE)
    if p_comb["ALLcos"] is not None:
        s1 = p_comb["ALLcos"]
        axLL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    axLL.set(xlim=(0, 0.7), xlabel=r"$x_{B}$", ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    ylo,yhi = auto_ylim([p_comb["ALLn0"], p_comb["ALLcos"]])
    axLL.set_ylim(ylo, yhi)
    axLL.grid(True, linestyle="--", alpha=0.6)
    harmLL = axLL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=0'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=1'),
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axLL.add_artist(harmLL)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{prefix}_CombinedOnly.pdf")
    plt.savefig(out_path)
    print(f"Saved combined-only figure: {out_path}")

# ----------- main ----------------------------------------------------
def main():
    if len(sys.argv) != 6:
        print("Usage: python plot_enpi_from_texts.py Su22.txt Fa22.txt Sp23.txt Combined.txt <Prefix>")
        print("Example: python plot_enpi_from_texts.py su22.txt fa22.txt sp23.txt combined.txt enpiHightGE")
        sys.exit(1)

    su22_path, fa22_path, sp23_path, comb_path, prefix = sys.argv[1:6]

    # Parse all four files
    su22 = parse_asym_file(su22_path)
    fa22 = parse_asym_file(fa22_path)
    sp23 = parse_asym_file(sp23_path)
    comb = parse_asym_file(comb_path)

    # Build per-period dicts
    p_su22 = build_period_dict(su22, prefix)
    p_fa22 = build_period_dict(fa22, prefix)
    p_sp23 = build_period_dict(sp23, prefix)
    p_comb = build_period_dict(comb,  prefix)

    # Output directory (keep your existing convention)
    out_dir = os.path.join("output", "enpi+")
    os.makedirs(out_dir, exist_ok=True)

    # Make the two figures
    plot_all_periods(p_su22, p_fa22, p_sp23, prefix, out_dir)
    plot_combined_only(p_comb, prefix, out_dir)

if __name__ == "__main__":
    main()