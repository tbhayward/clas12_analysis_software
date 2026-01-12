#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ep -> en pi+ asymmetries versus -t' in several xB bins, for three run periods
and for the combined (inverse-variance weighted) file.

Usage:
  python plot_enpi_from_texts.py Su22.txt Fa22.txt Sp23.txt Combined.txt [--rad signed_delta_summary.txt] [--migration migration_delta_pairs.txt]

Outputs (for each available xB bin):
  output/enpi+/rgc_enpi_<BinTag>_AllPeriods.pdf
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly.pdf

If --rad and/or --migration are provided, also:
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly_withRad.pdf
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly_withMig.pdf
  output/enpi+/rgc_enpi_<BinTag>_CombinedOnly_withRadAndMig.pdf

Overlay (Combined-only, 1x3):
  output/enpi+/rgc_enpi_xBOverlay_CombinedOnly.pdf
  (and if --rad and/or --migration)
    output/enpi+/rgc_enpi_xBOverlay_CombinedOnly_withRad.pdf
    output/enpi+/rgc_enpi_xBOverlay_CombinedOnly_withMig.pdf
    output/enpi+/rgc_enpi_xBOverlay_CombinedOnly_withRadAndMig.pdf

Notes:
- Fit-result text files must contain blocks like:
    enpiLowxBGEchi2FitsALUsinphi = {{mean_tprime, value, error}, ...};
  x is plotted as -t'.

- With --rad:
    * central values are shifted by the signed Delta from the summary (y <- y + Delta)
    * statistical errors remain unchanged
    * the radiative systematic is shown as a rectangle spanning the full fixed bin,
      from y=0 to y=sigma_rad at each -t' bin

- With --migration:
    * no central-value shift is applied
    * a bin-migration systematic is shown using the "_delta" lists in the provided file:
        enpi<BinTag>chi2Fits<Thing>_delta = {{mean_tprime, delta_value}, ...};
      the size of the migration systematic is taken as abs(delta_value)
    * this migration systematic is also drawn as a full-bin rectangle from y=0 to y=sigma_mig

- If both --rad and --migration are provided:
    * both components are drawn
    * the total systematic used for final banding is the quadrature sum:
        sigma_tot = sqrt(sigma_rad^2 + sigma_mig^2)
      and is drawn as an additional full-bin rectangle

- IMPORTANT (systematic bars & binning):
  The systematic rectangles span full fixed bins, not centered on data means.
  Edges (given to you in -t, we use their absolute values on the +(-t') axis):
    * For enpiGE (integrated), 0.10-wide bins with edges:
        0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25
    * For xB sub-bins, 0.20-wide bins with edges:
        0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25

- TSA systematics use only F_{UL}^{sin phi}; sin2 phi systematics are ignored.
- DSA systematics use only F_{LL}; F_{LL}^{cos phi} systematics are ignored.
- On the 1x3 overlay, systematic bands (if any) are taken only from the lowest xB bin (enpiLowxBGE).
"""

import sys
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch


# ---------------------------------------------------------------------
# Styling / knobs
# ---------------------------------------------------------------------
COLORS = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green",
    "Combined": "black",
}

# Colors for the four xB slices on the 1x3 overlay
XB_COLORS = {
    "enpiLowxBGE":     "tab:blue",
    "enpiMidLowxBGE":  "tab:orange",
    "enpiMidHighxBGE": "tab:green",
    "enpiHighxBGE":    "tab:red",
}

# Marker shapes per xB slice (overlay)
XB_MARKERS = {
    "enpiLowxBGE":     "o",
    "enpiMidLowxBGE":  "s",
    "enpiMidHighxBGE": "^",
    "enpiHighxBGE":    "D",
}

MARKER = "o"
CAPSIZE = 3
MS = 5.0  # marker size

# Axes ranges/labels (x is -t')
XLIM_T = (0.0, 1.30)
X_LABEL = r"$-t'\ (\mathrm{GeV}^{2})$"

YLIM_LU = (-0.3, 0.3)   # BSA
YLIM_UL = (-0.3, 0.3)   # TSA
YLIM_LL = (-0.8, 0.8)   # DSA
YLIM_UU = (-0.4, 0.4)   # UU

# Human-readable labels for xB bins (for legends & titles)
XB_BINS = {
    "enpiLowxBGE":     r"$0.10 < x_{B} < 0.25$",
    "enpiMidLowxBGE":  r"$0.25 < x_{B} < 0.35$",
    "enpiMidHighxBGE": r"$0.35 < x_{B} < 0.45$",
    "enpiHighxBGE":    r"$0.45 < x_{B} < 0.60$",
    "enpiGE":          r"$0.10 < x_{B} < 0.60$"
}

# Order to render canvases
BIN_ORDER = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE", "enpiGE"]

# Top padding
TOP_PAD_PER_BIN = 0.95   # for 2x2 canvases with a suptitle
TOP_PAD_OVERLAY = 0.94   # for 1x3 overlay (no visible title)

# Output prefix constant
OUT_PREFIX = "enpi"

# Fixed bin-edge sets for systematic rectangles (on +(-t') axis)
EDGES_SYS_GE = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25], dtype=float)
EDGES_SYS_XB = np.array([0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25], dtype=float)


# ---------------------------------------------------------------------
# Parsing helpers (fit-result .txt files: NAME = {{x,y,dy},...};)
# ---------------------------------------------------------------------
_assign_re = re.compile(r'([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};', re.DOTALL)
_brace_entry_re = re.compile(r'\{([^{}]+)\}')

def parse_asym_file(path):
    """Parse NAME = {{mean,val,err}, ...}; blocks into dict[name] -> np.array(N,3)."""
    with open(path, "r") as f:
        text = f.read()
    out = {}
    for m in _assign_re.finditer(text):
        name = m.group(1)
        content = m.group(2)
        triples = []
        for t in _brace_entry_re.findall(content):
            parts = [p.strip() for p in t.split(",")]
            if len(parts) != 3:
                continue
            try:
                triples.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
            #endif
        #endfor
        if triples:
            out[name] = np.array(triples, dtype=float)
        #endif
    #endfor
    return out


def get_series(dct, key, negate_x=True, sort_x=True):
    """
    Return dict(x,y,yerr) if key exists with finite values & positive errors; else None.
    Negates x if negate_x=True (to convert t' -> -t' for plotting). Optionally sort by x.
    """
    if key not in dct:
        return None
    #endif
    arr = np.array(dct[key], dtype=float)
    if arr.size == 0:
        return None
    #endif
    x_raw, y, e = arr[:, 0], arr[:, 1], arr[:, 2]
    mask = np.isfinite(x_raw) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    if not np.any(mask):
        return None
    #endif
    x = -x_raw[mask] if negate_x else x_raw[mask]
    y = y[mask]
    e = e[mask]
    if sort_x and x.size > 1:
        order = np.argsort(x)
        x, y, e = x[order], y[order], e[order]
    #endif
    return {"x": x, "y": y, "yerr": e}


def build_period_dict(parsed, bin_prefix):
    """
    Collect the series for one period given a bin_prefix like 'enpiLowxBGE' or 'enpiGE'.
    We look for keys of the form f"{bin_prefix}chi2Fits<suffix>".
    """
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}"
    return {
        # BSA
        "ALUsin":  get_series(parsed, k("ALUsinphi")),
        # TSA (n=1 open, n=2 filled)
        "AULsin":  get_series(parsed, k("AULsinphi")),
        "AULsin2": get_series(parsed, k("AULsin2phi")),
        # DSA (n=0 open, n=1 (cos phi) filled)
        "ALLn0":   get_series(parsed, k("ALL")),
        "ALLcos":  get_series(parsed, k("ALLcosphi")),
        # UU (n=1 open, n=2 filled)
        "UUcos":   get_series(parsed, k("AUUcosphi")),
        "UUcos2":  get_series(parsed, k("AUUcos2phi")),
    }


def detect_available_bins(*dicts):
    """
    Inspect provided parsed dictionaries and return the subset of BIN_ORDER whose
    expected keys exist in ANY of the dicts. We probe ALUsinphi to decide presence.
    """
    available = []
    for bin_tag in BIN_ORDER:
        key_probe = f"{bin_tag}chi2FitsALUsinphi"
        present = any((d is not None and key_probe in d) for d in dicts)
        if present:
            available.append(bin_tag)
        #endif
    #endfor
    return available


# ---------------------------------------------------------------------
# Migration systematic parser (pairs: NAME = {{x,delta},...};)
# ---------------------------------------------------------------------
def parse_pair_file(path):
    """Parse NAME = {{x,val}, ...}; blocks into dict[name] -> np.array(N,2)."""
    with open(path, "r") as f:
        text = f.read()
    out = {}
    for m in _assign_re.finditer(text):
        name = m.group(1)
        content = m.group(2)
        pairs = []
        for t in _brace_entry_re.findall(content):
            parts = [p.strip() for p in t.split(",")]
            if len(parts) != 2:
                continue
            try:
                pairs.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
            #endif
        #endfor
        if pairs:
            out[name] = np.array(pairs, dtype=float)
        #endif
    #endfor
    return out


def get_sigma_from_pairs(dct, key, negate_x=True, sort_x=True):
    """
    For a NAME = {{x,delta},...} list:
      return dict(x, sigma) where sigma = abs(delta)
    """
    if key not in dct:
        return None
    #endif
    arr = np.array(dct[key], dtype=float)
    if arr.size == 0:
        return None
    #endif
    x_raw, delta = arr[:, 0], arr[:, 1]
    mask = np.isfinite(x_raw) & np.isfinite(delta)
    if not np.any(mask):
        return None
    #endif
    x = -x_raw[mask] if negate_x else x_raw[mask]
    sigma = np.abs(delta[mask])
    if sort_x and x.size > 1:
        order = np.argsort(x)
        x, sigma = x[order], sigma[order]
    #endif
    return {"x": x, "sigma": sigma}


def build_migration_dict(parsed_pairs, bin_prefix):
    """
    Collect migration sigmas for one bin_prefix from the *_delta lists.
    Only the second entry in each pair is used, as abs(delta_value).
    """
    k = lambda suffix: f"{bin_prefix}chi2Fits{suffix}"
    return {
        "ALUsin":  get_sigma_from_pairs(parsed_pairs, k("ALUsinphi_delta")),
        "AULsin":  get_sigma_from_pairs(parsed_pairs, k("AULsinphi_delta")),
        "AULsin2": get_sigma_from_pairs(parsed_pairs, k("AULsin2phi_delta")),
        "ALLn0":   get_sigma_from_pairs(parsed_pairs, k("ALL_delta")),
        "ALLcos":  get_sigma_from_pairs(parsed_pairs, k("ALLcosphi_delta")),
    }


# ---------------------------------------------------------------------
# Radiative-correction summary parser (signed Delta + sigma_rad)
# ---------------------------------------------------------------------
# Map series by distinctive LaTeX substring present in the summary file.
_SERIES_KEYWORDS = {
    "ALUsin":   r"F_{LU}^{\sin\phi}/F_{UU}",
    "AULsin":   r"F_{UL}^{\sin\phi}/F_{UU}",
    "AULsin2":  r"F_{UL}^{\sin2\phi}/F_{UU}",
    "ALLn0":    r"F_{LL}/F_{UU}",
    "ALLcos":   r"F_{LL}^{\cos\phi}/F_{UU}",
}

def parse_rad_summary(rad_path):
    """
    Parse the signed-Delta summary text into:
      rad[bin_tag][series_key] = dict(x=array, delta=array, sigma=array)

    Expected structure (flexible):
      Bin: enpiLowxBGE ...
      Series: F_{LU}^{sin phi}/F_{UU}
      (optional header line)
      <x> <Delta> <sigma>
      ...
    """
    rad = {}
    cur_bin = None
    with open(rad_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("Bin:"):
            parts = line.split()
            if len(parts) >= 2:
                cur_bin = parts[1]
                rad.setdefault(cur_bin, {})
            #endif
            i += 1
            continue
        #endif

        if line.startswith("Series:"):
            ser_line = line[len("Series:"):].strip()
            cur_series_key = None
            for key, tag in _SERIES_KEYWORDS.items():
                if tag in ser_line:
                    cur_series_key = key
                    break
                #endif
            #endfor

            # Move to next line; if it looks like a header, skip it.
            i += 1
            if i < len(lines):
                maybe_hdr = lines[i].strip()
                if ("-t" in maybe_hdr) or ("Delta" in maybe_hdr) or ("sigma" in maybe_hdr):
                    i += 1
                #endif
            #endif

            xs, deltas, sigmas = [], [], []
            while i < len(lines):
                row = lines[i].strip()
                if (not row) or row.startswith("Series:") or row.startswith("Bin:") or (set(row) <= set("=-")):
                    break
                #endif
                parts = row.split()
                try:
                    xval = float(parts[0])
                    dval = float(parts[1])
                    sval = float(parts[2])
                except Exception:
                    break
                #endif
                xs.append(xval)
                deltas.append(dval)
                sigmas.append(sval)
                i += 1
            #endfor

            if cur_bin is not None and cur_series_key is not None and xs:
                rad[cur_bin][cur_series_key] = {
                    "x": np.array(xs, dtype=float),
                    "delta": np.array(deltas, dtype=float),
                    "sigma": np.array(sigmas, dtype=float),
                }
            #endif
            continue
        #endif

        i += 1
    #endfor

    return rad


# ---------------------------------------------------------------------
# Titles (2x2 canvases only)
# ---------------------------------------------------------------------
def make_title(xb_label):
    return r"$ep \rightarrow en\pi^{+}$ - " + xb_label


# ---------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------
def _suffix(with_rad, with_mig):
    if with_rad and with_mig:
        return "_withRadAndMig"
    #endif
    if with_rad:
        return "_withRad"
    #endif
    if with_mig:
        return "_withMig"
    #endif
    return ""


def _legend_run_period(ax, labels):
    handles = [Line2D([0], [0], marker=MARKER, color=COLORS[L], linestyle='', label=L) for L in labels]
    leg = ax.legend(handles=handles, title="Run Period", frameon=True, edgecolor="black",
                    loc="lower right", fontsize=11, title_fontsize=12)
    leg.get_frame().set_alpha(0.9)


def _legend_harmonic(ax, labels=("n=1", "n=2"), where="lower left"):
    open_marker = Line2D([0], [0], marker=MARKER, mfc='none', mec='black', linestyle='', label=labels[0])
    filled_marker = Line2D([0], [0], marker=MARKER, color='black', linestyle='', label=labels[1])
    leg = ax.legend(handles=[open_marker, filled_marker], title="Harmonic",
                    frameon=True, edgecolor="black", loc=where,
                    bbox_to_anchor=(0.02, 0.02), fontsize=11, title_fontsize=12)
    ax.add_artist(leg)


def _sys_edges_for_bin(bin_tag):
    return EDGES_SYS_GE if bin_tag == "enpiGE" else EDGES_SYS_XB


def _draw_sys_bands_fullbin(ax, sigma, edges, *, facecolor="0.7", alpha=0.35, hatch=None, edgecolor=None, zorder=1):
    """
    Draw systematic rectangles that span the entire fixed x-bin:
      x in [edges[i], edges[i+1]]
      y in [0, sigma[i]]

    The number of rectangles drawn is min(len(sigma), len(edges)-1).
    """
    if sigma is None or edges is None:
        return
    #endif
    sigma = np.asarray(sigma, dtype=float)
    edges = np.asarray(edges, dtype=float)
    nbins = max(0, len(edges) - 1)
    if nbins == 0 or sigma.size == 0:
        return
    #endif

    n = min(nbins, sigma.size)
    if sigma.size != nbins:
        print(f"[INFO] sys sigma count ({sigma.size}) != #bins ({nbins}); drawing first {n} bins.")
    #endif

    for i in range(n):
        sv = sigma[i]
        if not np.isfinite(sv):
            continue
        #endif
        x0 = edges[i]
        x1 = edges[i + 1]
        if (not np.isfinite(x0)) or (not np.isfinite(x1)) or (x1 <= x0):
            continue
        #endif
        rect = Rectangle(
            (x0, 0.0), x1 - x0, sv,
            facecolor=facecolor,
            edgecolor=edgecolor,
            alpha=alpha,
            hatch=hatch,
            zorder=zorder
        )
        ax.add_patch(rect)
    #endfor


def _apply_rad_shift_if_available(series, rad_entry):
    """
    Apply y <- y + Delta using index pairing (no nearest-neighbor matching).
    Returns (x, y_shift, yerr, sigma_rad_vector_or_none).
    """
    if series is None:
        return None
    #endif
    x, y, yerr = series["x"], series["y"], series["yerr"]

    if rad_entry is None:
        return x, y, yerr, None
    #endif

    xr = np.asarray(rad_entry["x"], dtype=float)
    dr = np.asarray(rad_entry["delta"], dtype=float)
    sr = np.asarray(rad_entry["sigma"], dtype=float)

    n = min(len(x), len(xr), len(dr), len(sr))
    if n == 0:
        return x, y, yerr, None
    #endif
    if len(x) != len(xr):
        print(f"[WARN] rad length ({len(xr)}) != data length ({len(x)}); pairing by index up to {n}.")
    #endif

    y_shift = y.copy()
    y_shift[:n] = y_shift[:n] + dr[:n]

    sigma_out = np.full_like(y, np.nan, dtype=float)
    sigma_out[:n] = sr[:n]

    return x, y_shift, yerr, sigma_out


def _extract_sigma_mig(mig_entry):
    if mig_entry is None:
        return None
    #endif
    return np.asarray(mig_entry["sigma"], dtype=float)


def _combine_sigma_quadrature(sigma_rad, sigma_mig):
    """
    Combine by index: sigma_tot[i] = sqrt(sigma_rad[i]^2 + sigma_mig[i]^2)
    If one component missing at i, sigma_tot takes the other.
    Returns array of length max(len(rad), len(mig)) with NaNs where neither exists.
    """
    sr = None if sigma_rad is None else np.asarray(sigma_rad, dtype=float)
    sm = None if sigma_mig is None else np.asarray(sigma_mig, dtype=float)

    n_rad = 0 if sr is None else sr.size
    n_mig = 0 if sm is None else sm.size
    nmax = max(n_rad, n_mig)

    if nmax == 0:
        return None
    #endif

    out = np.full(nmax, np.nan, dtype=float)
    for i in range(nmax):
        has_r = (sr is not None) and (i < n_rad) and np.isfinite(sr[i])
        has_m = (sm is not None) and (i < n_mig) and np.isfinite(sm[i])
        if has_r and has_m:
            out[i] = np.sqrt(sr[i] * sr[i] + sm[i] * sm[i])
        elif has_r:
            out[i] = sr[i]
        elif has_m:
            out[i] = sm[i]
        #endif
    #endfor

    return out


def _legend_systematics(ax, with_rad, with_mig):
    """
    Add a small legend indicating the systematic components drawn.
    Convention used here:
      - migration: hatched outline
      - radiative: hatched outline (different hatch)
      - total (quadrature): filled gray rectangle
    """
    handles = []
    if with_rad and with_mig:
        handles.append(Patch(facecolor="0.75", edgecolor=None, alpha=0.35, label="sys (total)"))
    #endif
    if with_rad:
        handles.append(Patch(facecolor="none", edgecolor="0.35", hatch="////", label="sys (rad)"))
    #endif
    if with_mig:
        handles.append(Patch(facecolor="none", edgecolor="0.55", hatch="\\\\\\\\", label="sys (migration)"))
    #endif

    if not handles:
        return
    #endif

    leg = ax.legend(handles=handles, frameon=True, edgecolor="black", fontsize=10,
                    loc="upper left", bbox_to_anchor=(0.02, 0.98))
    leg.get_frame().set_alpha(0.9)
    ax.add_artist(leg)


def _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label,
                     *, bin_tag_for_width,
                     with_rad=False, rad_for_bin=None,
                     with_mig=False, mig_for_bin=None):
    """
    Draw the four panels for a set of labeled period dicts.

    Systematics:
      - radiative: uses rad_for_bin[series_key]["sigma"] and also shifts y by Delta
      - migration: uses mig_for_bin[series_key]["sigma"] = abs(delta_value) and does NOT shift y
      - if both: total = quadrature(rad, migration) is also drawn

    TSA: only AULsin systematic bands are used (ignore sin2).
    DSA: only ALLn0 systematic bands are used (ignore ALLcos).
    UU: no systematic bands are drawn here (consistent with the five ratios of interest).
    """
    edges = _sys_edges_for_bin(bin_tag_for_width)

    def R(key):
        if rad_for_bin is None:
            return None
        #endif
        return rad_for_bin.get(key)
    #endif

    def M(key):
        if mig_for_bin is None:
            return None
        #endif
        return mig_for_bin.get(key)
    #endif

    # ---------------- BSA (ALUsin) ----------------
    sigma_rad = None
    if with_rad and (R("ALUsin") is not None):
        sigma_rad = np.asarray(R("ALUsin")["sigma"], dtype=float)
    #endif
    sigma_mig = _extract_sigma_mig(M("ALUsin")) if with_mig else None
    sigma_tot = _combine_sigma_quadrature(sigma_rad, sigma_mig) if (with_rad and with_mig) else None

    if with_rad and with_mig and (sigma_tot is not None):
        _draw_sys_bands_fullbin(axLU, np.abs(sigma_tot), edges, facecolor="0.75", alpha=0.35, hatch=None, edgecolor=None, zorder=1)
    #endif
    if with_rad and (sigma_rad is not None):
        _draw_sys_bands_fullbin(axLU, np.abs(sigma_rad), edges, facecolor="none", alpha=1.0, hatch="////", edgecolor="0.35", zorder=2)
    #endif
    if with_mig and (sigma_mig is not None):
        _draw_sys_bands_fullbin(axLU, np.abs(sigma_mig), edges, facecolor="none", alpha=1.0, hatch="\\\\\\\\", edgecolor="0.55", zorder=3)
    #endif

    for lab, pdata in pdata_by_label.items():
        s = pdata["ALUsin"]
        if s is None:
            continue
        #endif
        if with_rad and (R("ALUsin") is not None):
            x, y, ye, _ = _apply_rad_shift_if_available(s, R("ALUsin"))
        else:
            x, y, ye = s["x"], s["y"], s["yerr"]
        #endif
        axLU.errorbar(x, y, yerr=ye, fmt=MARKER, color=COLORS[lab], ecolor=COLORS[lab],
                      capsize=CAPSIZE, label=lab, markersize=MS, linestyle="None")
    #endfor
    axLU.set(xlim=XLIM_T, ylim=YLIM_LU, xlabel=X_LABEL, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
    axLU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLU.grid(True, linestyle="--", alpha=0.6)
    _legend_run_period(axLU, list(pdata_by_label.keys()))
    _legend_systematics(axLU, with_rad=with_rad, with_mig=with_mig)

    # ---------------- TSA (AULsin open, AULsin2 filled) ----------------
    sigma_rad = None
    if with_rad and (R("AULsin") is not None):
        sigma_rad = np.asarray(R("AULsin")["sigma"], dtype=float)
    #endif
    sigma_mig = _extract_sigma_mig(M("AULsin")) if with_mig else None
    sigma_tot = _combine_sigma_quadrature(sigma_rad, sigma_mig) if (with_rad and with_mig) else None

    if with_rad and with_mig and (sigma_tot is not None):
        _draw_sys_bands_fullbin(axUL, np.abs(sigma_tot), edges, facecolor="0.75", alpha=0.35, hatch=None, edgecolor=None, zorder=1)
    #endif
    if with_rad and (sigma_rad is not None):
        _draw_sys_bands_fullbin(axUL, np.abs(sigma_rad), edges, facecolor="none", alpha=1.0, hatch="////", edgecolor="0.35", zorder=2)
    #endif
    if with_mig and (sigma_mig is not None):
        _draw_sys_bands_fullbin(axUL, np.abs(sigma_mig), edges, facecolor="none", alpha=1.0, hatch="\\\\\\\\", edgecolor="0.55", zorder=3)
    #endif

    for lab, pdata in pdata_by_label.items():
        s1 = pdata["AULsin"]
        if s1 is not None:
            if with_rad and (R("AULsin") is not None):
                x1, y1, ye1, _ = _apply_rad_shift_if_available(s1, R("AULsin"))
            else:
                x1, y1, ye1 = s1["x"], s1["y"], s1["yerr"]
            #endif
            axUL.errorbar(x1, y1, yerr=ye1, fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif

        s2 = pdata["AULsin2"]
        if s2 is not None:
            # No TSA systematic bars from sin2, but we still plot points; optionally apply rad shift if available.
            if with_rad and (R("AULsin2") is not None):
                x2, y2, ye2, _ = _apply_rad_shift_if_available(s2, R("AULsin2"))
            else:
                x2, y2, ye2 = s2["x"], s2["y"], s2["yerr"]
            #endif
            axUL.errorbar(x2, y2, yerr=ye2, fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif
    #endfor
    axUL.set(xlim=XLIM_T, ylim=YLIM_UL, xlabel=X_LABEL, ylabel=r"$F_{UL}^{\sin n\phi}/F_{UU}$")
    axUL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUL.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axUL, labels=("n=1", "n=2"))
    _legend_run_period(axUL, list(pdata_by_label.keys()))
    _legend_systematics(axUL, with_rad=with_rad, with_mig=with_mig)

    # ---------------- DSA (ALLn0 open, ALLcos filled) ----------------
    sigma_rad = None
    if with_rad and (R("ALLn0") is not None):
        sigma_rad = np.asarray(R("ALLn0")["sigma"], dtype=float)
    #endif
    sigma_mig = _extract_sigma_mig(M("ALLn0")) if with_mig else None
    sigma_tot = _combine_sigma_quadrature(sigma_rad, sigma_mig) if (with_rad and with_mig) else None

    if with_rad and with_mig and (sigma_tot is not None):
        _draw_sys_bands_fullbin(axLL, np.abs(sigma_tot), edges, facecolor="0.75", alpha=0.35, hatch=None, edgecolor=None, zorder=1)
    #endif
    if with_rad and (sigma_rad is not None):
        _draw_sys_bands_fullbin(axLL, np.abs(sigma_rad), edges, facecolor="none", alpha=1.0, hatch="////", edgecolor="0.35", zorder=2)
    #endif
    if with_mig and (sigma_mig is not None):
        _draw_sys_bands_fullbin(axLL, np.abs(sigma_mig), edges, facecolor="none", alpha=1.0, hatch="\\\\\\\\", edgecolor="0.55", zorder=3)
    #endif

    for lab, pdata in pdata_by_label.items():
        s0 = pdata["ALLn0"]
        if s0 is not None:
            if with_rad and (R("ALLn0") is not None):
                x0, y0, ye0, _ = _apply_rad_shift_if_available(s0, R("ALLn0"))
            else:
                x0, y0, ye0 = s0["x"], s0["y"], s0["yerr"]
            #endif
            axLL.errorbar(x0, y0, yerr=ye0, fmt="o", mfc="none", mec=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif

        s1 = pdata["ALLcos"]
        if s1 is not None:
            # No DSA systematic bars from cos; we still plot points; optionally apply rad shift if available.
            if with_rad and (R("ALLcos") is not None):
                x1, y1, ye1, _ = _apply_rad_shift_if_available(s1, R("ALLcos"))
            else:
                x1, y1, ye1 = s1["x"], s1["y"], s1["yerr"]
            #endif
            axLL.errorbar(x1, y1, yerr=ye1, fmt="o", color=COLORS[lab], ecolor=COLORS[lab],
                          capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif
    #endfor
    axLL.set(xlim=XLIM_T, ylim=YLIM_LL, xlabel=X_LABEL, ylabel=r"$F_{LL}^{\cos n\phi}/F_{UU}$")
    axLL.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axLL.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axLL, labels=("n=0", "n=1"))
    _legend_run_period(axLL, list(pdata_by_label.keys()))
    _legend_systematics(axLL, with_rad=with_rad, with_mig=with_mig)

    # ---------------- UU (no sys bands) ----------------
    for lab, pdata in pdata_by_label.items():
        s1 = pdata["UUcos"]
        if s1 is not None:
            axUU.errorbar(s1["x"], s1["y"], yerr=s1["yerr"], fmt="o", mfc="none",
                          mec=COLORS[lab], ecolor=COLORS[lab], capsize=CAPSIZE,
                          markersize=MS, linestyle="None")
        #endif
        s2 = pdata["UUcos2"]
        if s2 is not None:
            axUU.errorbar(s2["x"], s2["y"], yerr=s2["yerr"], fmt="o", color=COLORS[lab],
                          ecolor=COLORS[lab], capsize=CAPSIZE, markersize=MS, linestyle="None")
        #endif
    #endfor
    axUU.set(xlim=XLIM_T, ylim=YLIM_UU, xlabel=X_LABEL, ylabel=r"$F_{UU}^{\cos n\phi}/F_{UU}$")
    axUU.axhline(0, color="black", linestyle="--", linewidth=1.2)
    axUU.grid(True, linestyle="--", alpha=0.6)
    _legend_harmonic(axUU, labels=("n=1", "n=2"))
    _legend_run_period(axUU, list(pdata_by_label.keys()))


# ---------------------------------------------------------------------
# Per-bin canvases (2x2)
# ---------------------------------------------------------------------
def plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, out_dir):
    plt.figure(figsize=(12, 9))  # 2x2
    axLU = plt.subplot(2, 2, 1)  # BSA
    axUL = plt.subplot(2, 2, 2)  # TSA
    axLL = plt.subplot(2, 2, 3)  # DSA
    axUU = plt.subplot(2, 2, 4)  # UU

    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(xb_label), fontsize=16, y=0.97)

    pdata_by_label = {"Su22": p_su22, "Fa22": p_fa22, "Sp23": p_sp23}
    _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label,
                     bin_tag_for_width=bin_tag,
                     with_rad=False, rad_for_bin=None,
                     with_mig=False, mig_for_bin=None)

    plt.tight_layout(rect=[0, 0, 1, TOP_PAD_PER_BIN])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_{bin_tag}_AllPeriods.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved all-periods figure: {out_path}")


def plot_combined_only_for_bin(p_comb, bin_tag, out_dir, *, with_rad, rad_bin, with_mig, mig_bin):
    plt.figure(figsize=(12, 9))  # 2x2
    axLU = plt.subplot(2, 2, 1)
    axUL = plt.subplot(2, 2, 2)
    axLL = plt.subplot(2, 2, 3)
    axUU = plt.subplot(2, 2, 4)

    xb_label = XB_BINS.get(bin_tag, bin_tag)
    plt.suptitle(make_title(xb_label), fontsize=16, y=0.97)

    pdata_by_label = {"Combined": p_comb}
    _plot_panel_sets(axLU, axUL, axLL, axUU, pdata_by_label,
                     bin_tag_for_width=bin_tag,
                     with_rad=with_rad, rad_for_bin=rad_bin,
                     with_mig=with_mig, mig_for_bin=mig_bin)

    plt.tight_layout(rect=[0, 0, 1, TOP_PAD_PER_BIN])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_{bin_tag}_CombinedOnly{_suffix(with_rad, with_mig)}.pdf")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved combined-only figure: {out_path}")


# ---------------------------------------------------------------------
# xB overlay 1x3 canvas (Combined only) - NO suptitle
# ---------------------------------------------------------------------
def plot_combined_xb_overlay_1x3(comb_parsed, out_dir, bins_to_use,
                                 *, with_rad, rad_all,
                                 with_mig, mig_all):
    """
    1x3 canvas (Combined only) overlaying the four xB slices:
      Left  : F_LU^{sin phi}/F_UU  (ALUsinphi)
      Middle: F_UL^{sin phi}/F_UU  (AULsinphi)
      Right : F_UL^{sin2phi}/F_UU  (AULsin2phi)

    Systematic rectangles on overlay:
      - taken ONLY from enpiLowxBGE
      - widths use the xB 0.20-wide fixed edges
      - TSA systematics: use only AULsin
      - if both rad and migration are available: also draw total quadrature band
    """
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharex=True)
    axL, axM, axR = axes

    low_bin_tag = "enpiLowxBGE"
    edges_overlay = EDGES_SYS_XB

    def _rad_entry(bin_tag, key):
        if (not with_rad) or (rad_all is None):
            return None
        #endif
        return rad_all.get(bin_tag, {}).get(key)
    #endif

    def _mig_entry(bin_tag, key):
        if (not with_mig) or (mig_all is None):
            return None
        #endif
        return mig_all.get(bin_tag, {}).get(key)
    #endif

    def _draw_sys_for_overlay(ax, series_key, draw_sys):
        if not draw_sys:
            return
        #endif

        # Components always taken from the lowest xB bin.
        r = _rad_entry(low_bin_tag, series_key)
        m = _mig_entry(low_bin_tag, series_key)

        sigma_rad = None
        if r is not None:
            sigma_rad = np.asarray(r["sigma"], dtype=float)
        #endif
        sigma_mig = None
        if m is not None:
            sigma_mig = np.asarray(m["sigma"], dtype=float)
        #endif
        sigma_tot = _combine_sigma_quadrature(sigma_rad, sigma_mig) if (sigma_rad is not None and sigma_mig is not None) else None

        if (sigma_tot is not None) and (sigma_rad is not None) and (sigma_mig is not None):
            _draw_sys_bands_fullbin(ax, np.abs(sigma_tot), edges_overlay, facecolor="0.75", alpha=0.35, hatch=None, edgecolor=None, zorder=1)
        #endif
        if sigma_rad is not None:
            _draw_sys_bands_fullbin(ax, np.abs(sigma_rad), edges_overlay, facecolor="none", alpha=1.0, hatch="////", edgecolor="0.35", zorder=2)
        #endif
        if sigma_mig is not None:
            _draw_sys_bands_fullbin(ax, np.abs(sigma_mig), edges_overlay, facecolor="none", alpha=1.0, hatch="\\\\\\\\", edgecolor="0.55", zorder=3)
        #endif

        _legend_systematics(ax, with_rad=(sigma_rad is not None), with_mig=(sigma_mig is not None))
    #endif

    def _draw_component(ax, series_key, ylabel, ylim, *, draw_sys):
        # Systematics (if requested for this component)
        _draw_sys_for_overlay(ax, series_key, draw_sys)

        # Points for each xB slice
        for bin_tag in bins_to_use:
            pdata = build_period_dict(comb_parsed, bin_tag)
            s = pdata[series_key]
            if s is None:
                continue
            #endif

            clr = XB_COLORS.get(bin_tag, "black")
            mkr = XB_MARKERS.get(bin_tag, "o")

            # Apply rad shift if present for this bin/series
            rad_entry = _rad_entry(bin_tag, series_key)
            if with_rad and (rad_entry is not None):
                x, y, ye, _ = _apply_rad_shift_if_available(s, rad_entry)
            else:
                x, y, ye = s["x"], s["y"], s["yerr"]
            #endif

            ax.errorbar(x, y, yerr=ye, fmt=mkr, color=clr, ecolor=clr,
                        capsize=CAPSIZE, label=XB_BINS.get(bin_tag, bin_tag),
                        markersize=MS, linestyle="None")
        #endfor

        ax.set(xlim=XLIM_T, ylim=ylim, xlabel=X_LABEL, ylabel=ylabel)
        ax.axhline(0, color="black", linestyle="--", linewidth=1.2)
        ax.grid(True, linestyle="--", alpha=0.6)

        leg = ax.legend(frameon=True, edgecolor="black", fontsize=11, loc="best")
        leg.get_frame().set_alpha(0.9)
    #endif

    # Left: ALUsin (draw sys)
    _draw_component(axL, "ALUsin", r"$F_{LU}^{\sin\phi}/F_{UU}$", YLIM_LU, draw_sys=True)
    # Middle: AULsin (draw sys; TSA uses only sin phi)
    _draw_component(axM, "AULsin", r"$F_{UL}^{\sin\phi}/F_{UU}$", YLIM_UL, draw_sys=True)
    # Right: AULsin2 (NO sys; by instruction we ignore sin2 phi systematics)
    _draw_component(axR, "AULsin2", r"$F_{UL}^{\sin2\phi}/F_{UU}$", YLIM_UL, draw_sys=False)

    fig.tight_layout(rect=[0, 0, 1, TOP_PAD_OVERLAY])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rgc_{OUT_PREFIX}_xBOverlay_CombinedOnly{_suffix(with_rad, with_mig)}.pdf")
    plt.savefig(out_path)
    plt.close(fig)
    print(f"Saved xB-overlay (1x3) figure: {out_path}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Plot enpi+ asymmetries vs -t' with optional radiative shifts and bin-migration systematics.")
    ap.add_argument("su22", type=str, help="Su22 fit-results .txt")
    ap.add_argument("fa22", type=str, help="Fa22 fit-results .txt")
    ap.add_argument("sp23", type=str, help="Sp23 fit-results .txt")
    ap.add_argument("combined", type=str, help="Combined fit-results .txt")
    ap.add_argument("--rad", type=str, default=None, help="Signed Delta summary .txt (applies y <- y + Delta and draws sigma_rad bars)")
    ap.add_argument("--migration", type=str, default=None, help="Migration delta-pairs .txt (uses *_delta lists; sigma_mig = abs(delta_value))")
    args = ap.parse_args()

    # Parse fit-result files
    su22 = parse_asym_file(args.su22)
    fa22 = parse_asym_file(args.fa22)
    sp23 = parse_asym_file(args.sp23)
    comb = parse_asym_file(args.combined)

    # Which xB-bin canvases can we make?
    available_bins = detect_available_bins(su22, fa22, sp23, comb)
    if not available_bins:
        print("[ERROR] No recognizable xB-bin sections found (e.g. enpiLowxBGEchi2FitsALUsinphi).")
        sys.exit(2)
    #endif

    # Optional radiative-correction summary (FAIL FAST if requested and missing)
    rad_all = None
    if args.rad is not None:
        if not os.path.isfile(args.rad):
            print(f"[ERROR] --rad file not found: {args.rad}")
            sys.exit(2)
        #endif
        rad_all = parse_rad_summary(args.rad)
    #endif

    # Optional migration systematics (FAIL FAST if requested and missing)
    mig_all = None
    mig_pairs = None
    if args.migration is not None:
        if not os.path.isfile(args.migration):
            print(f"[ERROR] --migration file not found: {args.migration}")
            sys.exit(2)
        #endif
        mig_pairs = parse_pair_file(args.migration)

        # Build per-bin migration dicts once, so plotting just does lookups.
        mig_all = {}
        for bin_tag in available_bins:
            mig_all[bin_tag] = build_migration_dict(mig_pairs, bin_tag)
        #endfor
    #endif

    out_dir = os.path.join("output", "enpi+")
    os.makedirs(out_dir, exist_ok=True)

    # Per-bin canvases
    for bin_tag in available_bins:
        p_su22 = build_period_dict(su22, bin_tag)
        p_fa22 = build_period_dict(fa22, bin_tag)
        p_sp23 = build_period_dict(sp23, bin_tag)
        p_comb = build_period_dict(comb, bin_tag)

        # Always make the nominal "AllPeriods" canvas (3 runs)
        plot_all_periods_for_bin(p_su22, p_fa22, p_sp23, bin_tag, out_dir)

        # Always make the nominal "CombinedOnly" canvas
        plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                   with_rad=False, rad_bin=None,
                                   with_mig=False, mig_bin=None)

        # If migration requested: CombinedOnly_withMig
        if mig_all is not None:
            plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                       with_rad=False, rad_bin=None,
                                       with_mig=True, mig_bin=mig_all.get(bin_tag, {}))
        #endif

        # If rad requested: CombinedOnly_withRad (and withRadAndMig if mig also present)
        if rad_all is not None:
            rad_bin = rad_all.get(bin_tag, {})
            plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                       with_rad=True, rad_bin=rad_bin,
                                       with_mig=False, mig_bin=None)

            if mig_all is not None:
                plot_combined_only_for_bin(p_comb, bin_tag, out_dir,
                                           with_rad=True, rad_bin=rad_bin,
                                           with_mig=True, mig_bin=mig_all.get(bin_tag, {}))
            #endif
        #endif
    #endfor

    # Combined-only xB overlay (1x3)
    xb_overlay_candidates = ["enpiLowxBGE", "enpiMidLowxBGE", "enpiMidHighxBGE", "enpiHighxBGE"]
    bins_to_use = [b for b in xb_overlay_candidates if b in available_bins]
    if bins_to_use:
        # Nominal
        plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use,
                                     with_rad=False, rad_all=None,
                                     with_mig=False, mig_all=None)

        # With migration only
        if mig_all is not None:
            plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use,
                                         with_rad=False, rad_all=None,
                                         with_mig=True, mig_all=mig_all)
        #endif

        # With rad only, and rad+mig
        if rad_all is not None:
            plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use,
                                         with_rad=True, rad_all=rad_all,
                                         with_mig=False, mig_all=None)

            if mig_all is not None:
                plot_combined_xb_overlay_1x3(comb, out_dir, bins_to_use,
                                             with_rad=True, rad_all=rad_all,
                                             with_mig=True, mig_all=mig_all)
            #endif
        #endif
    else:
        print("[WARN] No xB-slice data available for overlay canvas; skipping.")
    #endif


if __name__ == "__main__":
    main()
#endif