#!/usr/bin/env python3
# apply_mx2_fermi_migration_to_asymmetries.py
#
# Purpose:
#   Estimate a bin-migration systematic in reconstructed Mx2 due to Fermi motion.
#
# Concept:
#   - You provide multiple fit-results text files, each produced with a different reconstructed Mx2 window.
#   - You also provide an Mx2 migration/composition matrix from MC.
#   - We identify the "exclusive reconstructed bin" (baseline window, typically 0.81 < Mx2 < 1.00).
#   - We take the migration row for that bin: weights w_k describing how much each Mx2 source-bin contributes
#     to the reconstructed exclusive bin.
#   - For each asymmetry series and each -t' point, we compute:
#         A_mix = sum_k w_k * A_k
#     where A_k is the fitted asymmetry from the fit file corresponding to Mx2 bin k.
#   - The migration systematic is:
#         Delta = A_mix - A_base
#         Sys   = |Delta|
#
# Key update (Jan 2026):
#   - The script NO LONGER treats small differences in the per-bin tmean grids across Mx2 windows as fatal.
#   - The mixing is performed by INDEX (bin number) assuming bin boundaries match; the output x-axis (tmean)
#     is taken from the baseline (exclusive) file.
#   - If the tmean grids differ by more than --tmean-tol, the script emits WARNING(s), not FATAL.
#
# Additional update (Jan 2026):
#   - Adds "Mx2-bin overlay" plots: for each xB group, a 2x3 grid (5 SF panels + 1 blank legend panel)
#     overlaying ALL Mx2-bin fit curves on each SF subplot.
#
# Outputs (written under output/enpi+/mx2_fermi_migration_study/):
#   (1) Migrated/mixed fit file:
#         name = {{tmean, A_mix, stat_base}, ...};
#   (2) Diff/systematics file:
#         name        = {{tmean, |A_mix - A_base|}, ...};
#         name_ratio  = {{tmean, |A_mix - A_base|/|A_base|}, ...};
#         name_delta  = {{tmean,  A_mix - A_base }, ...};
#   (3) PDF plots:
#         - polarized_structure_functions_baseline_<xBgroup>.pdf
#         - polarized_structure_functions_migrated_<xBgroup>.pdf
#         - migration_systematics_<xBgroup>.pdf
#         - mx2_bin_overlays_<xBgroup>.pdf
#
# Notes:
#   - This script supports dropping the final Mx2 bin (e.g., if you have 7 bins in the CSV but only 6 fit files).
#     If you provide N-1 fit files for an N-bin CSV, you MUST pass --drop-last-bin.
#   - The mixing uses only ONE migration row: the row corresponding to the reconstructed exclusive Mx2 bin.
#   - The baseline statistical uncertainties are preserved (stat_base).
#
# ------------------------------------------------------------------------------

import argparse
import ast
import csv
import os
import re
import sys


def fatal(msg):
    sys.stderr.write("FATAL: " + msg + "\n")
    sys.exit(1)
#enddef


def warn(msg):
    sys.stderr.write("WARNING: " + msg + "\n")
#enddef


def ensure_output_dir(out_dir):
    try:
        os.makedirs(out_dir, exist_ok=True)
    except Exception as e:
        fatal("Failed creating output directory '{}': {}".format(out_dir, str(e)))
    #endif
#enddef


def ensure_is_basename(path_like, what):
    if os.path.basename(path_like) != path_like:
        fatal("{} must be a file name only (no directory). Got: '{}'".format(what, path_like))
    #endif
#enddef


def read_mx2_migration_csv(csv_path):
    """
    Reads an Mx2 migration/composition CSV with columns:
      bin, mx2min, mx2max, nev, frac0, frac1, ..., frac(N-1)

    Returns:
      fracs_by_row: list[list[float]]  shape (N, N)
      meta:         list[(bin, mx2min, mx2max, nev)]
    """
    if not os.path.isfile(csv_path):
        fatal("Migration CSV not found: " + csv_path)
    #endif

    rows = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None:
            fatal("Migration CSV is empty: " + csv_path)
        #endif

        for row_idx, row in enumerate(reader):
            if len(row) < 5:
                fatal("Migration CSV row {} has too few columns (need at least 5).".format(row_idx))
            #endif

            try:
                bin_num = int(float(row[0]))
                mx2min = float(row[1])
                mx2max = float(row[2])
                nev = int(float(row[3]))
            except Exception as e:
                fatal("Failed parsing meta columns on CSV row {}: {}".format(row_idx, str(e)))
            #endif

            frac_cols = row[4:]
            fracs = []
            for j in range(len(frac_cols)):
                try:
                    fracs.append(float(frac_cols[j]))
                except Exception as e:
                    fatal("Failed parsing fraction col j={} on CSV row {}: {}".format(j, row_idx, str(e)))
                #endif
            #endfor

            rows.append((bin_num, mx2min, mx2max, nev, fracs))
        #endfor
    #endwith

    if len(rows) == 0:
        fatal("Migration CSV has no data rows: " + csv_path)
    #endif

    # Determine N and validate square
    N = len(rows)
    fracs_by_row = []
    meta = []

    for i in range(N):
        (bin_num, mx2min, mx2max, nev, fracs) = rows[i]
        if len(fracs) != N:
            fatal(
                "Migration CSV is not square: row {} has {} fraction columns, but CSV has {} rows.".format(
                    i, len(fracs), N
                )
            )
        #endif
        fracs_by_row.append(fracs)
        meta.append((bin_num, mx2min, mx2max, nev))
    #endfor

    # Validate bin numbering (must match row index; fail fast)
    for i in range(N):
        if meta[i][0] != i:
            fatal("Migration CSV bin mismatch: row {} has bin {} (expected {}).".format(i, meta[i][0], i))
        #endif
    #endfor

    return fracs_by_row, meta
#enddef


def find_exclusive_bin_index(meta, mx2min_target, mx2max_target, tol):
    """
    Find bin index whose (mx2min, mx2max) matches target within tol.
    """
    candidates = []
    for i in range(len(meta)):
        mx2min = meta[i][1]
        mx2max = meta[i][2]
        if (abs(mx2min - mx2min_target) <= tol) and (abs(mx2max - mx2max_target) <= tol):
            candidates.append(i)
        #endif
    #endfor

    if len(candidates) == 0:
        msg = "Could not find exclusive Mx2 bin matching [{:.6f}, {:.6f}] within tol {:.3g}.\n".format(
            mx2min_target, mx2max_target, tol
        )
        msg += "Available bins from CSV:\n"
        for i in range(len(meta)):
            msg += "  bin {}: [{:.6f}, {:.6f}]  nev={}\n".format(i, meta[i][1], meta[i][2], meta[i][3])
        #endfor
        fatal(msg.rstrip("\n"))
    #endif

    if len(candidates) > 1:
        fatal("Multiple CSV bins match exclusive window; ambiguous matches: {}".format(candidates))
    #endif

    return candidates[0]
#enddef


def parse_fit_results_text(txt_path):
    if not os.path.isfile(txt_path):
        fatal("Fit-results text file not found: " + txt_path)
    #endif

    text = open(txt_path, "r").read()

    pattern = re.compile(r"([A-Za-z0-9_]+)\s*=\s*(\{\{.*?\}\})\s*;", re.DOTALL)
    matches = list(pattern.finditer(text))
    if len(matches) == 0:
        fatal("No Mathematica-style assignments found in: " + txt_path)
    #endif

    out = {}
    for m in matches:
        name = m.group(1).strip()
        rhs = m.group(2).strip()

        py_rhs = rhs.replace("{", "[").replace("}", "]")
        try:
            val = ast.literal_eval(py_rhs)
        except Exception as e:
            fatal("Failed to parse assignment for variable '{}' in file '{}': {}".format(name, txt_path, str(e)))
        #endif

        out[name] = val
    #endfor

    return out
#enddef


def xb_group_name_from_index(idx):
    if idx == 0:
        return "enpiLowxB"
    elif idx == 1:
        return "enpiMidLowxB"
    elif idx == 2:
        return "enpiMidHighxB"
    elif idx == 3:
        return "enpiHighxB"
    else:
        fatal("Invalid xB group index: {}".format(idx))
    #endif
#enddef


def build_required_varnames():
    xb_groups = ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]
    suffixes = [
        "GEchi2FitsALUsinphi",
        "GEchi2FitsAULsinphi",
        "GEchi2FitsAULsin2phi",
        "GEchi2FitsALL",
        "GEchi2FitsALLcosphi",
    ]
    varnames = []
    for g in xb_groups:
        for s in suffixes:
            varnames.append(g + s)
        #endfor
    #endfor
    return varnames
#enddef


def validate_fit_map(fit_map, required_varnames, expected_len, file_tag):
    missing = []
    for v in required_varnames:
        if v not in fit_map:
            missing.append(v)
        #endif
    #endfor

    if len(missing) > 0:
        msg = "Missing required variables in fit-results file '{}':\n".format(file_tag)
        for v in missing:
            msg += "  " + v + "\n"
        #endfor
        fatal(msg.rstrip("\n"))
    #endif

    for v in required_varnames:
        arr = fit_map[v]
        if (not isinstance(arr, list)) or (len(arr) != expected_len):
            fatal(
                "Variable '{}' in file '{}' must be a list of length {}; got type {} length {}.".format(
                    v, file_tag, expected_len, type(arr).__name__, len(arr) if isinstance(arr, list) else -1
                )
            )
        #endif

        for k in range(expected_len):
            tup = arr[k]
            if (not isinstance(tup, (list, tuple))) or (len(tup) != 3):
                fatal(
                    "Variable '{}' in file '{}' entry {} must be a triple (tmean, value, stat).".format(
                        v, file_tag, k
                    )
                )
            #endif
        #endfor
    #endfor
#enddef


def warn_tmean_grids_if_inconsistent(fit_maps, required_varnames, expected_len, tmean_tol):
    """
    Non-fatal check.

    Mixing is done by INDEX; output/plot x-grid uses baseline tmean grid.

    Behavior:
      - If abs(tmean_k - tmean_ref) > tmean_tol, emit WARNING.
      - No fatal exit.
    """
    if len(fit_maps) < 2:
        return
    #endif

    ref = fit_maps[0]

    worst_abs = 0.0
    worst_info = None
    n_warn = 0

    for v in required_varnames:
        for i in range(expected_len):
            t_ref = float(ref[v][i][0])
            for k in range(1, len(fit_maps)):
                t_k = float(fit_maps[k][v][i][0])
                dt = abs(t_k - t_ref)
                if dt > tmean_tol:
                    n_warn += 1
                    if n_warn <= 25:
                        warn(
                            "tmean grid differs (non-fatal) for '{}' idx {}: file0 {:.9f} vs file{} {:.9f} (abs diff {:.3g} > tol {:.3g}).".format(
                                v, i, t_ref, k, t_k, dt, tmean_tol
                            )
                        )
                    #endif
                #endif
                if dt > worst_abs:
                    worst_abs = dt
                    worst_info = (v, i, t_ref, k, t_k, dt)
                #endif
            #endfor
        #endfor
    #endfor

    if n_warn > 25:
        warn("tmean grid: emitted first 25 warnings; suppressed {} additional warnings.".format(n_warn - 25))
    #endif

    if worst_info is not None and worst_abs > tmean_tol:
        (v, i, t_ref, k, t_k, dt) = worst_info
        warn(
            "tmean grid summary: max abs diff {:.3g} occurred for '{}' idx {}: file0 {:.9f} vs file{} {:.9f}.".format(
                dt, v, i, t_ref, k, t_k
            )
        )
    #endif
#enddef


def compute_mx2_mixed_values(weights_row, fit_maps, excl_bin_idx, required_varnames, renormalize, tol):
    """
    weights_row: list[float] weights for each provided fit file (length M)
    fit_maps:    list[dict]  one per Mx2 bin file (length M)
    excl_bin_idx: int index of the exclusive/baseline file in fit_maps list

    Mixing is performed by INDEX. Output tmean is taken from the baseline file.
    """
    if len(weights_row) != len(fit_maps):
        fatal("Internal error: weights_row length {} != number of fit_maps {}.".format(len(weights_row), len(fit_maps)))
    #endif

    row_sum = 0.0
    for w in weights_row:
        row_sum += w
    #endfor

    if row_sum == 0.0:
        fatal("Exclusive-bin migration row sum is zero; cannot compute weighted mixture.")
    #endif

    if (abs(row_sum - 1.0) > tol) and (not renormalize):
        fatal(
            "Exclusive-bin migration row sum is {:.6f}, not within tol {:.6f} of 1.0.\n"
            "If you intend to renormalize the weights, rerun with --renormalize.".format(
                row_sum, tol
            )
        )
    #endif

    weights = list(weights_row)
    if renormalize:
        weights = [w / row_sum for w in weights_row]
    #endif

    migrated = {}
    diffs_mag = {}
    diffs_signed = {}

    for v in required_varnames:
        base_list = fit_maps[excl_bin_idx][v]

        out_triples = []
        out_diff_mag_pairs = []
        out_diff_signed_pairs = []

        for i in range(len(base_list)):
            tmean = float(base_list[i][0])     # always use baseline tmean grid
            aval_base = float(base_list[i][1])
            stat_base = float(base_list[i][2]) # preserve baseline stat

            aval_mix = 0.0
            for k in range(len(fit_maps)):
                a_k = float(fit_maps[k][v][i][1])  # index-paired mixing
                aval_mix += weights[k] * a_k
            #endfor

            delta = aval_mix - aval_base
            sysmag = abs(delta)

            out_triples.append([tmean, aval_mix, stat_base])
            out_diff_mag_pairs.append([tmean, sysmag])
            out_diff_signed_pairs.append([tmean, delta])
        #endfor

        migrated[v] = out_triples
        diffs_mag[v] = out_diff_mag_pairs
        diffs_signed[v] = out_diff_signed_pairs
    #endfor

    return migrated, diffs_mag, diffs_signed
#enddef


def format_mathematica_list_triple(lst):
    parts = []
    for trip in lst:
        parts.append("{" + "{:.9f}".format(trip[0]) + ", " + "{:.9f}".format(trip[1]) + ", " + "{:.9f}".format(trip[2]) + "}")
    #endfor
    return "{" + ", ".join(parts) + "}"
#enddef


def format_mathematica_list_pair(lst):
    parts = []
    for pair in lst:
        parts.append("{" + "{:.9f}".format(pair[0]) + ", " + "{:.9f}".format(pair[1]) + "}")
    #endfor
    return "{" + ", ".join(parts) + "}"
#enddef


def write_output_files(out_path, out_diff_path, migrated, diffs_mag, diffs_signed, fit_map_baseline, required_varnames):
    with open(out_path, "w") as f:
        for v in required_varnames:
            rhs = format_mathematica_list_triple(migrated[v])
            f.write(v + " = " + rhs + ";\n")
        #endfor
        f.write("\n")
    #endwith

    eps = 1.0e-15

    with open(out_diff_path, "w") as f:
        for v in required_varnames:
            if v not in diffs_mag:
                fatal("Internal error: diffs_mag missing '{}'.".format(v))
            #endif
            if v not in diffs_signed:
                fatal("Internal error: diffs_signed missing '{}'.".format(v))
            #endif
            if v not in fit_map_baseline:
                fatal("Internal error: baseline fit map missing '{}'.".format(v))
            #endif

            rhs_mag = format_mathematica_list_pair(diffs_mag[v])
            f.write(v + " = " + rhs_mag + ";\n")

            rhs_signed = format_mathematica_list_pair(diffs_signed[v])
            f.write(v + "_delta = " + rhs_signed + ";\n")

            ratio_list = []
            base_list = fit_map_baseline[v]
            for i in range(len(base_list)):
                tmean = float(diffs_mag[v][i][0])
                sysmag = float(diffs_mag[v][i][1])

                before_val = float(base_list[i][1])
                denom = abs(before_val)
                if denom < eps:
                    fatal(
                        "Cannot compute ratio for {} at index {}: abs(baseline value) is {:.17g} (too close to zero).".format(
                            v, i, denom
                        )
                    )
                #endif

                ratio = sysmag / denom
                ratio_list.append([tmean, ratio])
            #endfor

            rhs_ratio = format_mathematica_list_pair(ratio_list)
            f.write(v + "_ratio = " + rhs_ratio + ";\n")

            f.write("\n")
        #endfor
    #endwith
#enddef


def default_diff_basename(out_basename):
    base, ext = os.path.splitext(out_basename)
    if ext == "":
        return out_basename + "_diff"
    #endif
    return base + "_diff" + ext
#enddef


# =========================
# Plotting
# =========================

def import_plot_deps():
    try:
        import numpy as np
    except Exception as e:
        fatal("numpy import failed (needed for plots): {}".format(str(e)))
    #endif

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        fatal("matplotlib import failed (needed for plots): {}".format(str(e)))
    #endif

    return np, plt
#enddef


def to_series(triples, negate_x=True, sort=True, np=None):
    if np is None:
        fatal("Internal error: numpy not provided to to_series().")
    #endif

    arr = np.array(triples, dtype=float)
    x = -arr[:, 0] if negate_x else arr[:, 0]
    y = arr[:, 1]
    e = arr[:, 2]
    if sort:
        idx = np.argsort(x)
        x, y, e = x[idx], y[idx], e[idx]
    #endif
    return x, y, e
#enddef


def compute_bin_edges_from_centers(x, np):
    n = len(x)
    if n < 2:
        fatal("Need at least 2 points to infer bin edges.")
    #endif

    mid = 0.5 * (x[1:] + x[:-1])

    left = np.zeros(n, dtype=float)
    right = np.zeros(n, dtype=float)

    left[1:] = mid
    right[:-1] = mid

    left[0] = x[0] - 0.5 * (x[1] - x[0])
    right[-1] = x[-1] + 0.5 * (x[-1] - x[-2])

    widths = right - left

    if np.any(widths <= 0.0):
        fatal("Inferred non-positive bin width(s) from x centers: {}".format(widths))
    #endif

    return left, right, widths
#enddef


def extract_sys_series_explicit(fit_map, sys_name, x_ref, np):
    raw = fit_map[sys_name]
    if (not isinstance(raw, list)) or (len(raw) != len(x_ref)):
        fatal("Systematics series '{}' exists but has length {} (expected {}).".format(
            sys_name, len(raw) if isinstance(raw, list) else -1, len(x_ref)
        ))
    #endif

    tmeans = []
    sysv = []
    for i in range(len(raw)):
        ent = raw[i]
        if (not isinstance(ent, (list, tuple))) or (len(ent) not in [2, 3]):
            fatal("Systematics series '{}' entry {} must be [tmean, sys] or [tmean, *, sys].".format(sys_name, i))
        #endif
        tmean = float(ent[0])
        if len(ent) == 2:
            s = float(ent[1])
        else:
            s = float(ent[2])
        #endif
        tmeans.append(tmean)
        sysv.append(abs(s))
    #endfor

    arr = np.array(list(zip(tmeans, sysv)), dtype=float)
    x = -arr[:, 0]
    s = arr[:, 1]
    idx = np.argsort(x)
    x = x[idx]
    s = s[idx]

    if (len(x) == len(x_ref)) and (not np.allclose(x, x_ref, rtol=0.0, atol=1.0e-6)):
        warn("x mismatch between '{}' systematics and stat series (index-paired).".format(sys_name))
    #endif

    return s
#enddef


def get_sys_band(fit_map, varname, x_ref, np):
    preferred = varname + "Sys"
    if preferred in fit_map:
        return extract_sys_series_explicit(fit_map, preferred, x_ref, np)
    #endif

    candidates = [
        varname + "Syst",
        varname + "SYS",
        varname + "SYST",
        varname + "_sys",
        varname + "_syst",
        varname + "_SYS",
        varname + "_SYST",
    ]
    for c in candidates:
        if c in fit_map:
            return extract_sys_series_explicit(fit_map, c, x_ref, np)
        #endif
    #endfor

    return None
#enddef


def draw_sys_bars(ax, x_centers, sys_band, widths):
    if sys_band is None:
        return
    #endif
    ax.bar(x_centers,  sys_band, width=widths, bottom=0.0, color="0.7", alpha=0.25, linewidth=0.0, align="center", zorder=0)
    ax.bar(x_centers, -sys_band, width=widths, bottom=0.0, color="0.7", alpha=0.25, linewidth=0.0, align="center", zorder=0)
#enddef


def draw_points(ax, x, y, e, label, color, marker, capsize, ms):
    ax.errorbar(
        x, y, yerr=e,
        fmt=marker, color=color, ecolor=color,
        capsize=capsize, markersize=ms, linestyle="None",
        label=label, zorder=2
    )
#enddef


def xb_label(group_name):
    labels = {
        "enpiLowxB":     r"$0.10 < x_{B} < 0.25$",
        "enpiMidLowxB":  r"$0.25 < x_{B} < 0.35$",
        "enpiMidHighxB": r"$0.35 < x_{B} < 0.45$",
        "enpiHighxB":    r"$0.45 < x_{B} < 0.60$",
    }
    return labels.get(group_name, group_name)
#enddef


def save_polarized_structure_function_canvases(fit_map, out_dir, file_prefix):
    """
    Makes the 1x3 polarized-structure-function style canvases:
      - LU
      - UL (sinphi, sin2phi)
      - LL (0, cosphi)

    Uses sys rectangles if VarSys exists.
    """
    np, plt = import_plot_deps()

    xlim_t = (0.0, 1.30)
    x_label = r"$-t'\ (\mathrm{GeV}^{2})$"

    ylim_single = (-0.40, 0.40)
    ylim_double = (-0.80, 0.80)

    suffix_lu = "GEchi2FitsALUsinphi"
    suffix_ul1 = "GEchi2FitsAULsinphi"
    suffix_ul2 = "GEchi2FitsAULsin2phi"
    suffix_ll0 = "GEchi2FitsALL"
    suffix_ll1 = "GEchi2FitsALLcosphi"

    groups = ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]

    for g in groups:
        v_lu = g + suffix_lu
        v_ul1 = g + suffix_ul1
        v_ul2 = g + suffix_ul2
        v_ll0 = g + suffix_ll0
        v_ll1 = g + suffix_ll1

        for v in [v_lu, v_ul1, v_ul2, v_ll0, v_ll1]:
            if v not in fit_map:
                fatal("Cannot make plots: missing required input series '{}' in fit-results map.".format(v))
            #endif
        #endfor

        x_lu, y_lu, e_lu = to_series(fit_map[v_lu], np=np)
        x_ul1, y_ul1, e_ul1 = to_series(fit_map[v_ul1], np=np)
        x_ul2, y_ul2, e_ul2 = to_series(fit_map[v_ul2], np=np)
        x_ll0, y_ll0, e_ll0 = to_series(fit_map[v_ll0], np=np)
        x_ll1, y_ll1, e_ll1 = to_series(fit_map[v_ll1], np=np)

        if not np.allclose(x_lu, x_ul1, rtol=0.0, atol=1.0e-10):
            fatal("Plot x mismatch: {} vs {} (bin {})".format(v_lu, v_ul1, g))
        #endif
        if not np.allclose(x_lu, x_ul2, rtol=0.0, atol=1.0e-10):
            fatal("Plot x mismatch: {} vs {} (bin {})".format(v_lu, v_ul2, g))
        #endif
        if not np.allclose(x_lu, x_ll0, rtol=0.0, atol=1.0e-10):
            fatal("Plot x mismatch: {} vs {} (bin {})".format(v_lu, v_ll0, g))
        #endif
        if not np.allclose(x_lu, x_ll1, rtol=0.0, atol=1.0e-10):
            fatal("Plot x mismatch: {} vs {} (bin {})".format(v_lu, v_ll1, g))
        #endif

        left, right, widths = compute_bin_edges_from_centers(x_lu, np)

        sys_lu = get_sys_band(fit_map, v_lu, x_lu, np)
        sys_ul1 = get_sys_band(fit_map, v_ul1, x_lu, np)
        sys_ul2 = get_sys_band(fit_map, v_ul2, x_lu, np)
        sys_ll0 = get_sys_band(fit_map, v_ll0, x_lu, np)
        sys_ll1 = get_sys_band(fit_map, v_ll1, x_lu, np)

        # Convention:
        #   - UL panel uses only lower harmonic sys (sinphi)
        #   - LL panel uses only lower harmonic sys (LL)
        sys_band_ul_panel = sys_ul1
        sys_band_ll_panel = sys_ll0

        if (sys_lu is None) and (sys_band_ul_panel is None) and (sys_band_ll_panel is None):
            warn("No systematics series found for bin '{}'; sys rectangles will be skipped.".format(g))
        #endif

        fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4))

        ax_left = axes[0]
        ax_mid = axes[1]
        ax_right = axes[2]

        draw_sys_bars(ax_left, x_lu, sys_lu, widths)
        draw_points(
            ax_left, x_lu, y_lu, e_lu,
            label=r"$F_{LU}^{\sin\phi}/F_{UU}$",
            color="tab:blue", marker="o", capsize=3, ms=5.0
        )
        ax_left.set(xlim=xlim_t, ylim=ylim_single, xlabel=x_label, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
        ax_left.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
        ax_left.grid(True, linestyle="--", alpha=0.6)

        draw_sys_bars(ax_mid, x_ul1, sys_band_ul_panel, widths)
        draw_points(
            ax_mid, x_ul1, y_ul1, e_ul1,
            label=r"$F_{UL}^{\sin\phi}/F_{UU}$",
            color="tab:red", marker="s", capsize=3, ms=5.0
        )
        draw_points(
            ax_mid, x_ul2, y_ul2, e_ul2,
            label=r"$F_{UL}^{\sin2\phi}/F_{UU}$",
            color="tab:green", marker="^", capsize=3, ms=5.0
        )
        ax_mid.set(xlim=xlim_t, ylim=ylim_single, xlabel=x_label, ylabel=r"$F_{UL}^{\sin(n\phi)}/F_{UU}$")
        ax_mid.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
        ax_mid.grid(True, linestyle="--", alpha=0.6)
        leg_mid = ax_mid.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_mid.get_frame().set_alpha(0.9)

        draw_sys_bars(ax_right, x_ll0, sys_band_ll_panel, widths)
        draw_points(
            ax_right, x_ll0, y_ll0, e_ll0,
            label=r"$F_{LL}/F_{UU}$",
            color="tab:purple", marker="o", capsize=3, ms=5.0
        )
        draw_points(
            ax_right, x_ll1, y_ll1, e_ll1,
            label=r"$F_{LL}^{\cos\phi}/F_{UU}$",
            color="tab:orange", marker="s", capsize=3, ms=5.0
        )
        ax_right.set(xlim=xlim_t, ylim=ylim_double, xlabel=x_label, ylabel=r"$F_{LL}^{(0,\cos\phi)}/F_{UU}$")
        ax_right.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
        ax_right.grid(True, linestyle="--", alpha=0.6)
        leg_right = ax_right.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_right.get_frame().set_alpha(0.9)

        plt.suptitle(r"$ep \rightarrow en\pi^{+}$" + " - " + xb_label(g), fontsize=14, y=0.98)
        plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])

        out_pdf = os.path.join(out_dir, "{}_{}.pdf".format(file_prefix, g))
        plt.savefig(out_pdf)
        plt.close(fig)

        sys.stdout.write("Saved plot: {}\n".format(out_pdf))
    #endfor
#enddef


def inject_migration_sys_for_plotting(fit_map, diffs_mag, required_varnames):
    out = dict(fit_map)

    for v in required_varnames:
        if v not in diffs_mag:
            fatal("Internal error: diffs_mag missing '{}' needed for plotting sys injection.".format(v))
        #endif

        sys_pairs = []
        for ent in diffs_mag[v]:
            if (not isinstance(ent, (list, tuple))) or (len(ent) != 2):
                fatal("Internal error: diffs_mag['{}'] entry must be [tmean, sysmag].".format(v))
            #endif
            tmean = float(ent[0])
            sysmag = float(ent[1])
            sys_pairs.append([tmean, sysmag])
        #endfor

        out[v + "Sys"] = sys_pairs
    #endfor

    return out
#enddef


def save_migration_sys_canvases(fit_map, out_dir):
    np, plt = import_plot_deps()

    xlim_t = (0.0, 1.30)
    x_label = r"$-t'\ (\mathrm{GeV}^{2})$"

    suffix_lu = "GEchi2FitsALUsinphi"
    suffix_ul1 = "GEchi2FitsAULsinphi"
    suffix_ul2 = "GEchi2FitsAULsin2phi"
    suffix_ll0 = "GEchi2FitsALL"
    suffix_ll1 = "GEchi2FitsALLcosphi"

    groups = ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]

    for g in groups:
        v_lu = g + suffix_lu
        v_ul1 = g + suffix_ul1
        v_ul2 = g + suffix_ul2
        v_ll0 = g + suffix_ll0
        v_ll1 = g + suffix_ll1

        v_lu_sys = v_lu + "Sys"
        v_ul1_sys = v_ul1 + "Sys"
        v_ul2_sys = v_ul2 + "Sys"
        v_ll0_sys = v_ll0 + "Sys"
        v_ll1_sys = v_ll1 + "Sys"

        for v in [v_lu, v_ul1, v_ul2, v_ll0, v_ll1]:
            if v not in fit_map:
                fatal("Cannot make migration-sys plots: missing required input series '{}'.".format(v))
            #endif
        #endfor

        for v in [v_lu_sys, v_ul1_sys, v_ul2_sys, v_ll0_sys, v_ll1_sys]:
            if v not in fit_map:
                fatal("Cannot make migration-sys plots: missing injected sys series '{}'.".format(v))
            #endif
        #endfor

        x_ref, y_dummy, e_dummy = to_series(fit_map[v_lu], np=np)

        sys_lu = get_sys_band(fit_map, v_lu, x_ref, np)
        sys_ul1 = get_sys_band(fit_map, v_ul1, x_ref, np)
        sys_ul2 = get_sys_band(fit_map, v_ul2, x_ref, np)
        sys_ll0 = get_sys_band(fit_map, v_ll0, x_ref, np)
        sys_ll1 = get_sys_band(fit_map, v_ll1, x_ref, np)

        if (sys_lu is None) or (sys_ul1 is None) or (sys_ul2 is None) or (sys_ll0 is None) or (sys_ll1 is None):
            fatal("Internal error: expected all sys series to exist for '{}' but at least one is missing.".format(g))
        #endif

        fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4))

        ax_lu = axes[0]
        ax_ul = axes[1]
        ax_ll = axes[2]

        y_lu_max = float(np.max(sys_lu))
        if y_lu_max <= 0.0:
            y_lu_max = 1.0e-6
        #endif
        ax_lu.plot(x_ref, sys_lu, marker="o", linestyle="None", label=r"$|\Delta(F_{LU}^{\sin\phi}/F_{UU})|$")
        ax_lu.set(xlim=xlim_t, ylim=(0.0, 1.2 * y_lu_max), xlabel=x_label, ylabel="Migration systematic magnitude")
        ax_lu.grid(True, linestyle="--", alpha=0.6)
        leg_lu = ax_lu.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_lu.get_frame().set_alpha(0.9)

        y_ul_max = float(np.max(np.array([np.max(sys_ul1), np.max(sys_ul2)], dtype=float)))
        if y_ul_max <= 0.0:
            y_ul_max = 1.0e-6
        #endif
        ax_ul.plot(x_ref, sys_ul1, marker="o", linestyle="None", label=r"$|\Delta(F_{UL}^{\sin\phi}/F_{UU})|$")
        ax_ul.plot(x_ref, sys_ul2, marker="s", linestyle="None", label=r"$|\Delta(F_{UL}^{\sin2\phi}/F_{UU})|$")
        ax_ul.set(xlim=xlim_t, ylim=(0.0, 1.2 * y_ul_max), xlabel=x_label, ylabel="Migration systematic magnitude")
        ax_ul.grid(True, linestyle="--", alpha=0.6)
        leg_ul = ax_ul.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_ul.get_frame().set_alpha(0.9)

        y_ll_max = float(np.max(np.array([np.max(sys_ll0), np.max(sys_ll1)], dtype=float)))
        if y_ll_max <= 0.0:
            y_ll_max = 1.0e-6
        #endif
        ax_ll.plot(x_ref, sys_ll0, marker="o", linestyle="None", label=r"$|\Delta(F_{LL}/F_{UU})|$")
        ax_ll.plot(x_ref, sys_ll1, marker="s", linestyle="None", label=r"$|\Delta(F_{LL}^{\cos\phi}/F_{UU})|$")
        ax_ll.set(xlim=xlim_t, ylim=(0.0, 1.2 * y_ll_max), xlabel=x_label, ylabel="Migration systematic magnitude")
        ax_ll.grid(True, linestyle="--", alpha=0.6)
        leg_ll = ax_ll.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_ll.get_frame().set_alpha(0.9)

        plt.suptitle(r"$ep \rightarrow en\pi^{+}$" + " - " + xb_label(g) + " - Mx2 migration systematics", fontsize=14, y=0.98)
        plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])

        out_pdf = os.path.join(out_dir, "migration_systematics_{}.pdf".format(g))
        plt.savefig(out_pdf)
        plt.close(fig)

        sys.stdout.write("Saved plot: {}\n".format(out_pdf))
    #endfor
#enddef


def save_mx2_bin_overlay_canvases(fit_maps, meta, effective_N, out_dir):
    """
    For each xB group, make a 2x3 canvas:
      - 5 subplots: each polarized structure function (as extracted from each Mx2-window file)
      - 1 blank subplot: used for legend

    Each SF subplot overlays all Mx2-bin curves (one curve per input fit file).
    """
    np, plt = import_plot_deps()

    xlim_t = (0.0, 1.30)
    x_label = r"$-t'\ (\mathrm{GeV}^{2})$"

    ylim_single = (-0.40, 0.40)
    ylim_double = (-0.80, 0.80)

    suffix_lu  = "GEchi2FitsALUsinphi"
    suffix_ul1 = "GEchi2FitsAULsinphi"
    suffix_ul2 = "GEchi2FitsAULsin2phi"
    suffix_ll0 = "GEchi2FitsALL"
    suffix_ll1 = "GEchi2FitsALLcosphi"

    groups = ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]

    mx2_bin_labels = []
    for k in range(effective_N):
        mx2min = meta[k][1]
        mx2max = meta[k][2]
        mx2_bin_labels.append("bin {}: [{:.3f}, {:.3f}]".format(k, mx2min, mx2max))
    #endfor

    for g in groups:
        v_lu  = g + suffix_lu
        v_ul1 = g + suffix_ul1
        v_ul2 = g + suffix_ul2
        v_ll0 = g + suffix_ll0
        v_ll1 = g + suffix_ll1

        fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.5))

        ax_lu = axes[0, 0]
        ax_ul1 = axes[0, 1]
        ax_ul2 = axes[0, 2]
        ax_ll0 = axes[1, 0]
        ax_ll1 = axes[1, 1]
        ax_blank = axes[1, 2]

        legend_handles = []
        legend_labels = []

        def overlay_one(ax, varname, ylim, collect_legend):
            nonlocal legend_handles, legend_labels

            for k in range(effective_N):
                fm = fit_maps[k]
                if varname not in fm:
                    fatal("Overlay plots: missing '{}' in fit map index {}.".format(varname, k))
                #endif

                x, y, e = to_series(fm[varname], np=np)

                h = ax.errorbar(x, y, yerr=e, fmt="o", capsize=2, markersize=4.0, linestyle="None")

                if collect_legend:
                    legend_handles.append(h[0])
                    legend_labels.append(mx2_bin_labels[k])
                #endif
            #endfor

            ax.set(xlim=xlim_t, ylim=ylim, xlabel=x_label)
            ax.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
            ax.grid(True, linestyle="--", alpha=0.6)
        #enddef

        overlay_one(ax_lu, v_lu, ylim_single, collect_legend=True)
        ax_lu.set_ylabel(r"$F_{LU}^{\sin\phi}/F_{UU}$")

        overlay_one(ax_ul1, v_ul1, ylim_single, collect_legend=False)
        ax_ul1.set_ylabel(r"$F_{UL}^{\sin\phi}/F_{UU}$")

        overlay_one(ax_ul2, v_ul2, ylim_single, collect_legend=False)
        ax_ul2.set_ylabel(r"$F_{UL}^{\sin2\phi}/F_{UU}$")

        overlay_one(ax_ll0, v_ll0, ylim_double, collect_legend=False)
        ax_ll0.set_ylabel(r"$F_{LL}/F_{UU}$")

        overlay_one(ax_ll1, v_ll1, ylim_double, collect_legend=False)
        ax_ll1.set_ylabel(r"$F_{LL}^{\cos\phi}/F_{UU}$")

        ax_blank.axis("off")
        if len(legend_handles) > 0:
            ax_blank.legend(
                legend_handles,
                legend_labels,
                loc="center",
                frameon=True,
                edgecolor="black",
                fontsize=10
            )
        #endif

        plt.suptitle(r"$ep \rightarrow en\pi^{+}$" + " - " + xb_label(g) + " - Mx2-bin overlays", fontsize=14, y=0.98)
        plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

        out_pdf = os.path.join(out_dir, "mx2_bin_overlays_{}.pdf".format(g))
        plt.savefig(out_pdf)
        plt.close(fig)

        sys.stdout.write("Saved plot: {}\n".format(out_pdf))
    #endfor
#enddef


def main():
    ap = argparse.ArgumentParser(
        description="Apply an Mx2 Fermi-motion migration mixture (from MC) to enpi* GEchi2Fits* asymmetry lists."
    )

    ap.add_argument("migration_csv", help="Mx2 migration CSV (square N x N with meta columns).")
    ap.add_argument(
        "fit_txts",
        nargs="+",
        help="List of Mathematica-style fit-results text files, one per Mx2 bin (in CSV bin order)."
    )
    ap.add_argument("fit_out_name", help="Output file name only (no directory). Saved under output/enpi+/mx2_fermi_migration_study/")

    ap.add_argument("--out-diff", default=None, help="Diff output file name only. Default: <fit_out_name with _diff inserted>")

    # Exclusive bin selection:
    ap.add_argument("--exclusive-bin", type=int, default=None, help="Exclusive reconstructed Mx2 bin index (overrides --exclusive-mx2-min/max).")
    ap.add_argument("--exclusive-mx2-min", type=float, default=0.81, help="Exclusive reconstructed Mx2 min for auto-match.")
    ap.add_argument("--exclusive-mx2-max", type=float, default=1.00, help="Exclusive reconstructed Mx2 max for auto-match.")
    ap.add_argument("--mx2-match-tol", type=float, default=1.0e-6, help="Tolerance for matching exclusive mx2min/mx2max in CSV.")

    # Weight handling:
    ap.add_argument("--renormalize", action="store_true", help="Renormalize exclusive-bin weights to sum to 1.0 (after any dropping).")
    ap.add_argument("--tol", type=float, default=1.0e-3, help="Tolerance for weight-sum checks against 1.0 (strict mode).")

    # Dropping last bin (explicit):
    ap.add_argument("--drop-last-bin", action="store_true", help="Drop the final Mx2 bin from the mixture (requires providing N-1 fit files).")
    ap.add_argument("--dropped-frac-warn", type=float, default=0.02, help="Warn if dropped last-bin fraction exceeds this value.")

    # Validation:
    ap.add_argument("--tmean-tol", type=float, default=1.0e-8, help="Tolerance for tmean-grid warnings across Mx2 fit files (non-fatal).")

    args = ap.parse_args()

    out_dir = "output/enpi+/mx2_fermi_migration_study"
    ensure_output_dir(out_dir)

    ensure_is_basename(args.fit_out_name, "fit_out_name")
    out_diff_name = args.out_diff
    if out_diff_name is None:
        out_diff_name = default_diff_basename(args.fit_out_name)
    #endif
    ensure_is_basename(out_diff_name, "out-diff")

    out_path = os.path.join(out_dir, args.fit_out_name)
    out_diff_path = os.path.join(out_dir, out_diff_name)

    fracs_by_row, meta = read_mx2_migration_csv(args.migration_csv)
    N = len(fracs_by_row)

    # Determine how many fit files are provided, enforce explicit dropping if needed.
    M = len(args.fit_txts)

    if M == N:
        if args.drop_last_bin:
            fatal("You passed --drop-last-bin but provided N fit files (N={}). Provide N-1 fit files when dropping.".format(N))
        #endif
        effective_N = N
    elif M == (N - 1):
        if not args.drop_last_bin:
            fatal(
                "You provided N-1 fit files ({}), but the CSV has N={} bins.\n"
                "If you intend to ignore the final Mx2 bin, rerun with --drop-last-bin.".format(M, N)
            )
        #endif
        effective_N = N - 1
    else:
        fatal("Number of fit files ({}) must be N ({}) or N-1 ({}).".format(M, N, N - 1))
    #endif

    # Identify exclusive bin index in the CSV (full-N index)
    if args.exclusive_bin is not None:
        excl_csv_idx = args.exclusive_bin
        if excl_csv_idx < 0 or excl_csv_idx >= N:
            fatal("--exclusive-bin {} out of range for N={} bins.".format(excl_csv_idx, N))
        #endif
    else:
        excl_csv_idx = find_exclusive_bin_index(meta, args.exclusive_mx2_min, args.exclusive_mx2_max, args.mx2_match_tol)
    #endif

    # Map exclusive bin index to fit-file index (must exist within effective_N)
    if excl_csv_idx >= effective_N:
        fatal(
            "Exclusive bin index {} is not available in the provided fit files (effective_N={}).\n"
            "You likely dropped the last bin but also selected it as exclusive; adjust selection.".format(excl_csv_idx, effective_N)
        )
    #endif

    excl_fit_idx = excl_csv_idx

    # Extract the exclusive-bin weights row from CSV
    full_row = fracs_by_row[excl_csv_idx]

    if args.drop_last_bin:
        dropped = full_row[N - 1]
        if dropped > args.dropped_frac_warn:
            warn(
                "Dropping final Mx2 bin (index {}) but its fraction in exclusive-bin row is {:.6f} (> warn threshold {}).".format(
                    N - 1, dropped, args.dropped_frac_warn
                )
            )
        #endif
        weights_row = full_row[:effective_N]
    else:
        weights_row = full_row[:effective_N]
    #endif

    # Parse all fit files
    required = build_required_varnames()

    fit_maps = []
    for i in range(effective_N):
        path = args.fit_txts[i]
        fm = parse_fit_results_text(path)
        validate_fit_map(fm, required, expected_len=6, file_tag=path)
        fit_maps.append(fm)
    #endfor

    # Non-fatal tmean warnings (mixing is index-based; output uses baseline tmean grid)
    warn_tmean_grids_if_inconsistent(fit_maps, required, expected_len=6, tmean_tol=args.tmean_tol)

    # Compute mixture/migrated and diffs
    migrated, diffs_mag, diffs_signed = compute_mx2_mixed_values(
        weights_row=weights_row,
        fit_maps=fit_maps,
        excl_bin_idx=excl_fit_idx,
        required_varnames=required,
        renormalize=args.renormalize,
        tol=args.tol
    )

    # Write outputs
    baseline_fit_map = fit_maps[excl_fit_idx]

    write_output_files(
        out_path=out_path,
        out_diff_path=out_diff_path,
        migrated=migrated,
        diffs_mag=diffs_mag,
        diffs_signed=diffs_signed,
        fit_map_baseline=baseline_fit_map,
        required_varnames=required
    )

    sys.stdout.write("Wrote migrated fit file: {}\n".format(out_path))
    sys.stdout.write("Wrote difference file:   {}\n".format(out_diff_path))

    # Plotting:
    #   - baseline points + migration sys rectangles
    #   - migrated points + migration sys rectangles
    baseline_for_plots = inject_migration_sys_for_plotting(baseline_fit_map, diffs_mag, required)
    migrated_for_plots = inject_migration_sys_for_plotting(migrated, diffs_mag, required)

    save_polarized_structure_function_canvases(
        baseline_for_plots, out_dir, file_prefix="polarized_structure_functions_baseline"
    )
    save_polarized_structure_function_canvases(
        migrated_for_plots, out_dir, file_prefix="polarized_structure_functions_migrated"
    )

    save_migration_sys_canvases(baseline_for_plots, out_dir)

    # New: overlay all Mx2-bin fit curves on a 2x3 grid per xB slice
    save_mx2_bin_overlay_canvases(fit_maps, meta, effective_N, out_dir)

    # Diagnostic summary of the exclusive row weights
    sys.stdout.write("\nExclusive reconstructed Mx2 bin (CSV index): {}\n".format(excl_csv_idx))
    sys.stdout.write("Exclusive reconstructed Mx2 window from CSV: [{:.6f}, {:.6f}]\n".format(meta[excl_csv_idx][1], meta[excl_csv_idx][2]))
    sys.stdout.write("Using {} fit files (effective_N={}); drop_last_bin={}\n".format(effective_N, effective_N, args.drop_last_bin))
    sys.stdout.write("Weights provided (NOTE: if --renormalize is set, the actual weights used are rescaled internally):\n")
    for i in range(effective_N):
        sys.stdout.write("  w_raw[{}] = {:.6f}\n".format(i, weights_row[i]))
    #endfor
#enddef


if __name__ == "__main__":
    main()
#endif