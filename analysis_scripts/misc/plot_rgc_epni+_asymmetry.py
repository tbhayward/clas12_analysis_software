#!/usr/bin/env python3
import sys
import os
import re
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

COLORS = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green",
    "Combined": "black",
}
MARKER = "o"
CAPSIZE = 3
LABEL_FONTSIZE = 13

# ----------- parsing -------------------------------------------------
_assign_re = re.compile(r'([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};', re.DOTALL)
_triple_re = re.compile(r'\{([^{}]+)\}')

def parse_asym_file(path):
    """Parse NAME = {{mean,val,err}, ...}; blocks into dict[name] -> np.array(N,3)."""
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
                triples.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
        if triples:
            out[name] = np.array(triples, dtype=float)
    return out

def get_series(dct, key):
    """Return dict(x,y,yerr) if key exists and has finite values with positive errors; else None."""
    if key not in dct:
        return None
    arr = np.array(dct[key], dtype=float)
    if arr.size == 0:
        return None
    x, y, e = arr[:,0], arr[:,1], arr[:,2]
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    if not np.any(mask):
        return None
    return {"x": x[mask], "y": y[mask], "yerr": e[mask]}

def build_period_dict(parsed, prefix):
    """Collect the series for one period given a prefix like 'enpiHightGE'."""
    key = lambda suffix: f"{prefix}chi2Fits{suffix}"
    return {
        # BSA
        "ALUsin":  get_series(parsed, key("ALUsinphi")),
        # TSA (n=1 open, n=2 filled)
        "AULsin":  get_series(parsed, key("AULsinphi")),
        "AULsin2": get_series(parsed, key("AULsin2phi")),
        # DSA (n=0 open, cosφ filled)
        "ALLn0":   get_series(parsed, key("ALL")),
        "ALLcos":  get_series(parsed, key("ALLcosphi")),
        # UU (n=1 open, n=2 filled)
        "UUcos":   get_series(parsed, key("AUUcosphi")),
        "UUcos2":  get_series(parsed, key("AUUcos2phi")),
    }

# ----------- plotting helpers ---------------------------------------
def plot_all_periods(p_su22, p_fa22, p_sp23, prefix, out_dir):
    plt.figure(figsize=(12, 9))
    plt.suptitle(
        r"$ep \rightarrow en\pi^{+}$, $0.07<|t|<0.7$, $z>0.55$, $y<0.65$, $0.75<M_{x}^{2}<1.05\ \mathrm{GeV}^{2}$",
        fontsize=16, y=0.97
    )

    # Order requested: BSA TL, TSA TR, DSA BL, UU BR
    axLU  = plt.subplot(2,2,1)  # BSA
    axUL  = plt.subplot(2,2,2)  # TSA
    axLL  = plt.subplot(2,2,3)  # DSA
    axUU  = plt.subplot(2,2,4)  # UU

    # ---- BSA (ALU sinφ) ----
    for lab, pdata in (("Su22", p_su22), ("Fa22", p_fa22), ("Sp23", p_sp23)):
        s = pdata["ALUsin"]
        if s is not None:
            axLU.errorbar(s["x"], s["y"], yerr=s["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=lab)
    axLU.set(xlim=(0, 0.7), ylim=(-0.4, 0.4),
             xlabel=r"$x_{B}$", ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)
    legLU = axLU.legend(title="Run Period", frameon=True, edgecolor="black",
                        fontsize=11, title_fontsize=12)
    legLU.get_frame().set_alpha(0.9)

    # ---- TSA (AUL sinφ open, sin2φ filled) ----
    for lab, pdata in (("Su22", p_su22), ("Fa22", p_fa22), ("Sp23", p_sp23)):
        s1 = pdata["AULsin"]
        if s1 is not None:
            axUL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt=MARKER, mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1")
        s2 = pdata["AULsin2"]
        if s2 is not None:
            axUL.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=2")
    axUL.set(xlim=(0, 0.7), ylim=(-0.4, 0.4),
             xlabel=r"$x_{B}$", ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    harmUL = axUL.legend(
        handles=[
            Line2D([0],[0], marker=MARKER, mfc='none', mec='black', linestyle='', label='n=1'),
            Line2D([0],[0], marker=MARKER, color='black', linestyle='', label='n=2')
        ],
        title="Harmonic", frameon=True, edgecolor="black", loc="upper right",
        fontsize=11, title_fontsize=12
    )
    axUL.add_artist(harmUL)
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

    # ---- DSA (ALL n=0 open, cosφ filled) ----
    for lab, pdata in (("Su22", p_su22), ("Fa22", p_fa22), ("Sp23", p_sp23)):
        s0 = pdata["ALLn0"]
        if s0 is not None:
            axLL.errorbar(s0["x"], s0["y"], yerr=s0["yerr"],
                          fmt=MARKER, mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=0")
        s1 = pdata["ALLcos"]
        if s1 is not None:
            axLL.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1")
    axLL.set(xlim=(0, 0.7), ylim=(-0.6, 0.6),
             xlabel=r"$x_{B}$", ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
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

    # ---- UU (cos nφ) - n=1 open, n=2 filled ----
    for lab, pdata in (("Su22", p_su22), ("Fa22", p_fa22), ("Sp23", p_sp23)):
        s1 = pdata["UUcos"]
        s2 = pdata["UUcos2"]
        if s1 is not None:
            axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"],
                          fmt=MARKER, mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=1")
        if s2 is not None:
            axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"],
                          fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, label=f"{lab}, n=2")
    axUU.set(xlim=(0, 0.7), ylim=(-1.0, 1.0),
             xlabel=r"$x_{B}$", ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
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

    # Order requested: BSA TL, TSA TR, DSA BL, UU BR
    axLU  = plt.subplot(2,2,1)
    axUL  = plt.subplot(2,2,2)
    axLL  = plt.subplot(2,2,3)
    axUU  = plt.subplot(2,2,4)

    black = COLORS["Combined"]

    # BSA
    if p_comb["ALUsin"] is not None:
        s = p_comb["ALUsin"]
        axLU.errorbar(s["x"], s["y"], yerr=s["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    axLU.set(xlim=(0, 0.7), ylim=(-0.4, 0.4),
             xlabel=r"$x_{B}$", ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
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
    axUL.set(xlim=(0, 0.7), ylim=(-0.4, 0.4),
             xlabel=r"$x_{B}$", ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
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
    axLL.set(xlim=(0, 0.7), ylim=(-0.8, 0.8),
             xlabel=r"$x_{B}$", ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
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

    # UU
    if p_comb["UUcos"] is not None:
        axUU.errorbar(p_comb["UUcos"]["x"], p_comb["UUcos"]["y"], yerr=p_comb["UUcos"]["yerr"],
                      fmt=MARKER, mfc="none", mec=black, ecolor=black, capsize=CAPSIZE)
    if p_comb["UUcos2"] is not None:
        axUU.errorbar(p_comb["UUcos2"]["x"], p_comb["UUcos2"]["y"], yerr=p_comb["UUcos2"]["yerr"],
                      fmt=MARKER, color=black, ecolor=black, capsize=CAPSIZE)
    axUU.set(xlim=(0, 0.7), ylim=(-1.0, 1.0),
             xlabel=r"$x_{B}$", ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
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

    su22 = parse_asym_file(su22_path)
    fa22 = parse_asym_file(fa22_path)
    sp23 = parse_asym_file(sp23_path)
    comb = parse_asym_file(comb_path)

    p_su22 = build_period_dict(su22, prefix)
    p_fa22 = build_period_dict(fa22, prefix)
    p_sp23 = build_period_dict(sp23, prefix)
    p_comb = build_period_dict(comb,  prefix)

    out_dir = os.path.join("output", "enpi+")
    os.makedirs(out_dir, exist_ok=True)

    plot_all_periods(p_su22, p_fa22, p_sp23, prefix, out_dir)
    plot_combined_only(p_comb, prefix, out_dir)

if __name__ == "__main__":
    main()