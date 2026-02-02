#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_full_t_enpiplus_structure_function_ratios.py

Read one or more text files containing Mathematica-style list assignments like:

  tempGEchi2FitsALUsinphi = {{tprime_mean, value, err}, {...}, ...};

There may be other lists in the same file(s); we ignore anything except the five
target list names:

  tempGEchi2FitsALUsinphi
  tempGEchi2FitsAULsinphi
  tempGEchi2FitsAULsin2phi
  tempGEchi2FitsALL
  tempGEchi2FitsALLcosphi

We make a 2x3 canvas. Each target list gets its own subplot (one subplot left
blank). We plot (-t') on the x-axis (i.e., x = -tprime_mean), with y = asymmetry
value and yerr = uncertainty.

Axes:
  - Single-spin (first 3): y in [-0.6, 0.6]
  - Double-spin (last 2):  y in [-1.0, 1.0]
  - x for all:             x in [0.0, 10.0]

Labels:
  y axes:
    "F_{LU}^{sin#phi}/F_{UU}"
    "F_{UL}^{sin#phi}/F_{UU}"
    "F_{UL}^{sin2#phi}/F_{UU}"
    "F_{LL}/F_{UU}"
    "F_{LL}^{cos#phi}/F_{UU}"
  x axis: "-t'"

Output:
  - Creates output/full_t_enpi+/ if needed.
  - Saves the 2x3 canvas to:
      output/full_t_enpi+/structure_function_ratios_vs_tprime.png

Optional ROOT:
  If the final command-line argument ends with ".root", we also open it with uproot,
  read the "PhysicsEvents" tree, and histogram the branch "tprime".
  That histogram is saved to:
      output/full_t_enpi+/tprime_hist.png

Usage examples:
  python plot_full_t_enpiplus_structure_function_ratios.py fits.txt
  python plot_full_t_enpiplus_structure_function_ratios.py fits1.txt fits2.txt
  python plot_full_t_enpiplus_structure_function_ratios.py fits.txt some.root
  python plot_full_t_enpiplus_structure_function_ratios.py fits1.txt fits2.txt some.root
"""

import ast
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt


TARGET_VARS = [
    "tempGEchi2FitsALUsinphi",
    "tempGEchi2FitsAULsinphi",
    "tempGEchi2FitsAULsin2phi",
    "tempGEchi2FitsALL",
    "tempGEchi2FitsALLcosphi",
]

Y_LABELS = {
    "tempGEchi2FitsALUsinphi":   "F_{LU}^{sin#phi}/F_{UU}",
    "tempGEchi2FitsAULsinphi":   "F_{UL}^{sin#phi}/F_{UU}",
    "tempGEchi2FitsAULsin2phi":  "F_{UL}^{sin2#phi}/F_{UU}",
    "tempGEchi2FitsALL":         "F_{LL}/F_{UU}",
    "tempGEchi2FitsALLcosphi":   "F_{LL}^{cos#phi}/F_{UU}",
}

SSA_VARS = set([
    "tempGEchi2FitsALUsinphi",
    "tempGEchi2FitsAULsinphi",
    "tempGEchi2FitsAULsin2phi",
])

DSA_VARS = set([
    "tempGEchi2FitsALL",
    "tempGEchi2FitsALLcosphi",
])


def read_all_text(paths):
    chunks = []
    for p in paths:
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            chunks.append(f.read())
        #endfor
    return "\n".join(chunks)


def find_assignment_block(text, varname):
    """
    Return the raw Mathematica brace block string for:
        varname = {{...}, {...}, ...};

    We do a simple scan:
      - find 'varname'
      - find '='
      - find first '{'
      - then brace-match until the corresponding closing brace is found
    """
    # Locate "varname" followed by optional whitespace and '='
    m = re.search(r"\b" + re.escape(varname) + r"\b\s*=", text)
    if m is None:
        return None
    #endif

    i = m.end()  # position after '='
    n = len(text)

    # Skip whitespace
    while i < n and text[i].isspace():
        i += 1
    #endfor

    # Expect a '{'
    if i >= n or text[i] != "{":
        return None
    #endif

    start = i
    depth = 0
    while i < n:
        c = text[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                # Include this closing brace
                end = i + 1
                return text[start:end]
            #endif
        #endif
        i += 1
    #endfor

    return None


def mathematica_braces_to_python_list(brace_block):
    """
    Convert a Mathematica-style brace block to a Python-evaluable list string:
      {{a,b,c},{d,e,f}}  ->  [[a,b,c],[d,e,f]]

    Then parse via ast.literal_eval.
    """
    if brace_block is None:
        return None
    #endif

    s = brace_block.replace("{", "[").replace("}", "]")
    try:
        obj = ast.literal_eval(s)
    except Exception:
        return None
    #endif
    return obj


def extract_target_lists(text):
    """
    Return dict: varname -> list of triplets [[tprime_mean, value, err], ...]
    """
    out = {}
    for v in TARGET_VARS:
        block = find_assignment_block(text, v)
        parsed = mathematica_braces_to_python_list(block)
        if parsed is not None:
            out[v] = parsed
        #endif
    #endfor
    return out


def coerce_triplets(triplets, varname):
    """
    Ensure we have numeric arrays x=-tprime, y, yerr.
    Skip malformed rows quietly (but keep deterministic behavior).
    """
    xs = []
    ys = []
    es = []

    for row in triplets:
        if not isinstance(row, (list, tuple)) or len(row) < 3:
            continue
        #endif
        try:
            tprime_mean = float(row[0])
            val = float(row[1])
            err = float(row[2])
        except Exception:
            continue
        #endif

        xs.append(-tprime_mean)  # plot -t'
        ys.append(val)
        es.append(err)
    #endfor

    if len(xs) == 0:
        raise RuntimeError("No valid triplets found for " + varname)
    #endif

    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)
    e = np.array(es, dtype=float)

    # Sort by x increasing for nicer lines (if used later)
    order = np.argsort(x)
    return x[order], y[order], e[order]


def make_canvas(data_map, outdir):
    """
    2x3 canvas. Five plots + one blank.
    """
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=False)
    axes = axes.flatten()

    # Desired subplot order (left-to-right, top-to-bottom)
    plot_order = [
        "tempGEchi2FitsALUsinphi",
        "tempGEchi2FitsAULsinphi",
        "tempGEchi2FitsAULsin2phi",
        "tempGEchi2FitsALL",
        "tempGEchi2FitsALLcosphi",
    ]

    for i in range(6):
        ax = axes[i]

        if i >= len(plot_order):
            ax.axis("off")
            continue
        #endif

        var = plot_order[i]
        if var not in data_map:
            ax.axis("off")
            ax.text(0.5, 0.5, "Missing: " + var, ha="center", va="center", transform=ax.transAxes)
            continue
        #endif

        x, y, e = coerce_triplets(data_map[var], var)

        ax.errorbar(
            x, y, yerr=e,
            fmt="o",
            capsize=3,
            markersize=4,
            linewidth=1.0,
        )

        ax.set_xlim(0.0, 10.0)

        if var in SSA_VARS:
            ax.set_ylim(-0.6, 0.6)
        elif var in DSA_VARS:
            ax.set_ylim(-1.0, 1.0)
        else:
            # Should not happen for our target vars, but keep deterministic behavior.
            ax.set_ylim(-1.0, 1.0)
        #endif

        ax.set_xlabel("-t'")
        ax.set_ylabel(Y_LABELS.get(var, var))
        ax.grid(True, alpha=0.3)
    #endfor

    fig.tight_layout()

    outpath = os.path.join(outdir, "structure_function_ratios_vs_tprime.png")
    fig.savefig(outpath, dpi=200)
    plt.close(fig)

    return outpath


def make_tprime_hist(root_path, outdir):
    """
    Optional histogram of tprime from ROOT file PhysicsEvents tree.
    """
    try:
        import uproot
    except ImportError as e:
        raise RuntimeError("uproot is required for ROOT input but is not installed: " + str(e))
    #endif

    f = uproot.open(root_path)
    if "PhysicsEvents" not in f:
        raise RuntimeError("ROOT file does not contain tree 'PhysicsEvents': " + root_path)
    #endif

    tree = f["PhysicsEvents"]
    if "tprime" not in tree.keys():
        raise RuntimeError("Tree 'PhysicsEvents' does not contain branch 'tprime': " + root_path)
    #endif

    arr = tree["tprime"].array(library="np")
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]

    if arr.size == 0:
        raise RuntimeError("No finite entries in tprime for: " + root_path)
    #endif

    # Use a robust range so huge outliers do not dominate the plot.
    qlo, qhi = np.quantile(arr, [0.01, 0.99])
    if not np.isfinite(qlo) or not np.isfinite(qhi) or qhi <= qlo:
        qlo = float(np.min(arr))
        qhi = float(np.max(arr))
    #endif

    nbins = 200
    bins = np.linspace(qlo, qhi, nbins + 1)

    fig = plt.figure(figsize=(10, 6))
    plt.hist(arr, bins=bins, histtype="step")
    plt.xlabel("tprime")
    plt.ylabel("Counts")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    outpath = os.path.join(outdir, "tprime_hist.png")
    plt.savefig(outpath, dpi=200)
    plt.close(fig)

    return outpath


def main():
    if len(sys.argv) < 2:
        print("ERROR: Need at least one input text file.", file=sys.stderr)
        print("Usage: python plot_full_t_enpiplus_structure_function_ratios.py <fits1.txt> [fits2.txt ...] [optional.root]", file=sys.stderr)
        sys.exit(2)
    #endif

    args = sys.argv[1:]

    root_path = None
    if len(args) >= 2 and args[-1].lower().endswith(".root"):
        root_path = args[-1]
        text_paths = args[:-1]
    else:
        text_paths = args
    #endif

    for p in text_paths:
        if not os.path.isfile(p):
            print("ERROR: Text file not found: " + p, file=sys.stderr)
            sys.exit(2)
        #endif
    #endfor

    if root_path is not None and (not os.path.isfile(root_path)):
        print("ERROR: ROOT file not found: " + root_path, file=sys.stderr)
        sys.exit(2)
    #endif

    outdir = os.path.join("output", "full_t_enpi+")
    os.makedirs(outdir, exist_ok=True)

    text = read_all_text(text_paths)
    data_map = extract_target_lists(text)

    missing = [v for v in TARGET_VARS if v not in data_map]
    if len(missing) > 0:
        print("WARNING: Missing one or more target lists:", file=sys.stderr)
        for v in missing:
            print("  - " + v, file=sys.stderr)
        #endfor
        print("Will still produce plots; missing ones will be blank/marked.", file=sys.stderr)
    #endif

    canvas_path = make_canvas(data_map, outdir)
    print("Wrote: " + canvas_path)

    if root_path is not None:
        hist_path = make_tprime_hist(root_path, outdir)
        print("Wrote: " + hist_path)
    #endif


if __name__ == "__main__":
    main()
#endif