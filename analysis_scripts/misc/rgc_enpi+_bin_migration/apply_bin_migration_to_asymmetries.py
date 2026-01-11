#!/usr/bin/env python3
# apply_bin_migration_to_asymmetries.py
#
# Usage (tcsh one-liner example):
# python3 apply_bin_migration_to_asymmetries.py mc_bin_migration.csv fit_results.txt fit_results_migrated.txt
#
# Notes:
# - The 3rd argument is a file name only (no directory). Output files are written to:
#     output/rgc_enpi+_bin_migration_study/
# - The script produces:
#     (1) Migrated fit file:
#         name = {{tmean, migrated_value, stat_original}, ...};
#     (2) Diff/systematics file:
#         name        = {{tmean, abs(migrated-original)}, ...};               (systematic magnitude)
#         name_ratio  = {{tmean, abs(migrated-original)/abs(original)}, ...}; (relative magnitude)
#         name_delta  = {{tmean, migrated-original}, ...};                    (signed delta, diagnostics)
# - The script also produces 4 PDF canvases (one per xB bin) in the same output directory,
#   plotting the input/original asymmetries from the fit-results text file:
#     left:   F_LU^{sin(phi)} / F_UU
#     center: F_UL^{sin(phi)} / F_UU   and   F_UL^{sin(2phi)} / F_UU
#     right:  F_LL / F_UU              and   F_LL^{cos(phi)} / F_UU
#
# Additional plots (added without changing previous functionality):
# - The script also produces 4 PDF canvases (one per xB bin) showing ONLY the
#   migration systematic magnitudes vs -t':
#     UL:  |Delta AUL sin(phi)| and |Delta AUL sin(2phi)|
#     LL:  |Delta ALL| and |Delta ALL cos(phi)|
#   The y-axis of each subplot is [0, 1.2 * max_plotted_in_that_subplot].
#
# Plot requirements (as requested):
# - Y ranges (main asymmetry canvases):
#     single-spin (LU, UL): [-0.4, 0.4]
#     double-spin (LL):     [-0.8, 0.8]
# - No horizontal error bars (no xerr on points).
# - Systematic uncertainty drawn as per-bin rectangles from y=0 to +/-sys,
#   spanning the inferred (-t) bin width.
# - UL panel: sys rectangles from lower-order harmonic only (AUL sin(phi)).
# - LL panel: sys rectangles from lower-order harmonic only (ALL).
#
# Implementation detail:
# - Migration-derived systematic magnitudes are injected into the plotting map
#   under the deterministic name "<varname>Sys" (pairs: [tmean, sysmag]).
#   This guarantees the plotter can draw systematics even if the input text file
#   contains no sys series.

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

def read_migration_csv(csv_path):
    if not os.path.isfile(csv_path):
        fatal("Migration CSV not found: " + csv_path)
    #endif

    weights = []
    meta = []

    with open(csv_path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None:
            fatal("Migration CSV is empty: " + csv_path)
        #endif

        for row_idx, row in enumerate(reader):
            if len(row) != 30:
                fatal("Migration CSV row {} has {} columns, expected 30 (6 meta + 24 fractions).".format(row_idx, len(row)))
            #endif

            try:
                bin_num = int(float(row[0]))
                xbmin   = float(row[1])
                xbmax   = float(row[2])
                tmin_n  = float(row[3])
                tmax_n  = float(row[4])
                nev     = int(float(row[5]))
            except Exception as e:
                fatal("Failed parsing meta columns on CSV row {}: {}".format(row_idx, str(e)))
            #endif

            fracs = []
            for j in range(24):
                try:
                    fracs.append(float(row[6 + j]))
                except Exception as e:
                    fatal("Failed parsing fraction column j={} on CSV row {}: {}".format(j, row_idx, str(e)))
                #endif
            #endfor

            weights.append(fracs)
            meta.append((bin_num, xbmin, xbmax, tmin_n, tmax_n, nev))
        #endfor
    #endwith

    if len(weights) != 24:
        fatal("Migration CSV has {} data rows; expected 24 (bins 0-23).".format(len(weights)))
    #endif

    for i in range(24):
        if meta[i][0] != i:
            fatal("Migration CSV binNum mismatch: row {} has binNum {} (expected {}).".format(i, meta[i][0], i))
        #endif
    #endfor

    return weights, meta
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
        rhs  = m.group(2).strip()

        py_rhs = rhs.replace("{", "[").replace("}", "]")
        try:
            val = ast.literal_eval(py_rhs)
        except Exception as e:
            fatal("Failed to parse assignment for variable '{}': {}".format(name, str(e)))
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

def validate_fit_data(fit_map, required_varnames):
    missing = []
    for v in required_varnames:
        if v not in fit_map:
            missing.append(v)
        #endif
    #endfor

    if len(missing) > 0:
        msg = "Missing required variables in fit-results file:\n"
        for v in missing:
            msg += "  " + v + "\n"
        #endfor
        fatal(msg.rstrip("\n"))
    #endif

    for v in required_varnames:
        arr = fit_map[v]
        if (not isinstance(arr, list)) or (len(arr) != 6):
            fatal("Variable '{}' must be a list of length 6; got type {} length {}.".format(
                v, type(arr).__name__, len(arr) if isinstance(arr, list) else -1
            ))
        #endif

        for k in range(6):
            tup = arr[k]
            if (not isinstance(tup, (list, tuple))) or (len(tup) != 3):
                fatal("Variable '{}' entry {} must be (tmean, value, stat).".format(v, k))
            #endif
        #endfor
    #endfor
#enddef

def global_bin_to_group_and_tindex(b):
    if b < 0 or b > 23:
        fatal("Global bin out of range: {}".format(b))
    #endif
    return (b // 6), (b % 6)
#enddef

def compute_migrated_values(weights, fit_map, required_varnames, renormalize, tol):
    suffixes = [
        "GEchi2FitsALUsinphi",
        "GEchi2FitsAULsinphi",
        "GEchi2FitsAULsin2phi",
        "GEchi2FitsALL",
        "GEchi2FitsALLcosphi",
    ]

    migrated = {}
    diffs_mag = {}
    diffs_signed = {}

    for suffix in suffixes:
        for xb_group_idx in range(4):
            xb_name = xb_group_name_from_index(xb_group_idx)
            varname = xb_name + suffix

            out_list = []
            diff_mag_list = []
            diff_signed_list = []

            for t_idx in range(6):
                target_global_bin = xb_group_idx * 6 + t_idx

                row = weights[target_global_bin]
                row_sum = sum(row)

                if row_sum == 0.0:
                    fatal("Migration row sum is zero for target bin {} (cannot compute weighted average).".format(target_global_bin))
                #endif

                if (abs(row_sum - 1.0) > tol) and (not renormalize):
                    fatal(
                        "Migration row sum for target bin {} is {:.6f}, not within tol {:.6f} of 1.0.\n"
                        "If you intend to renormalize rows, rerun with --renormalize.".format(
                            target_global_bin, row_sum, tol
                        )
                    )
                #endif

                tmean_target = float(fit_map[varname][t_idx][0])
                aval_target  = float(fit_map[varname][t_idx][1])
                stat_target  = float(fit_map[varname][t_idx][2])

                migrated_val = 0.0
                for j in range(24):
                    w = row[j]
                    if w == 0.0:
                        continue
                    #endif

                    src_group_idx, src_t_idx = global_bin_to_group_and_tindex(j)
                    src_group_name = xb_group_name_from_index(src_group_idx)
                    src_varname = src_group_name + suffix

                    a_src = float(fit_map[src_varname][src_t_idx][1])
                    migrated_val += w * a_src
                #endfor

                if renormalize:
                    migrated_val = migrated_val / row_sum
                #endif

                delta_signed = migrated_val - aval_target
                delta_mag = abs(delta_signed)

                out_list.append([tmean_target, migrated_val, stat_target])
                diff_mag_list.append([tmean_target, delta_mag])
                diff_signed_list.append([tmean_target, delta_signed])
            #endfor

            migrated[varname] = out_list
            diffs_mag[varname] = diff_mag_list
            diffs_signed[varname] = diff_signed_list
        #endfor
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

def write_output_files(out_path, out_diff_path, migrated, diffs_mag, diffs_signed, fit_map):
    suffixes = [
        "GEchi2FitsALUsinphi",
        "GEchi2FitsAULsinphi",
        "GEchi2FitsAULsin2phi",
        "GEchi2FitsALL",
        "GEchi2FitsALLcosphi",
    ]
    xb_groups = ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]

    with open(out_path, "w") as f:
        for g in xb_groups:
            for s in suffixes:
                name = g + s
                if name not in migrated:
                    fatal("Internal error: migrated missing variable '{}'".format(name))
                #endif
                rhs = format_mathematica_list_triple(migrated[name])
                f.write(name + " = " + rhs + ";\n")
            #endfor
            f.write("\n")
        #endfor
    #endwith

    eps = 1.0e-15

    with open(out_diff_path, "w") as f:
        for g in xb_groups:
            for s in suffixes:
                name = g + s
                if name not in diffs_mag:
                    fatal("Internal error: diffs_mag missing variable '{}'".format(name))
                #endif
                if name not in diffs_signed:
                    fatal("Internal error: diffs_signed missing variable '{}'".format(name))
                #endif
                if name not in fit_map:
                    fatal("Internal error: fit_map missing variable '{}' needed for ratio.".format(name))
                #endif

                rhs_mag = format_mathematica_list_pair(diffs_mag[name])
                f.write(name + " = " + rhs_mag + ";\n")

                rhs_signed = format_mathematica_list_pair(diffs_signed[name])
                f.write(name + "_delta = " + rhs_signed + ";\n")

                ratio_list = []
                for k in range(6):
                    tmean = float(diffs_mag[name][k][0])
                    sysmag = float(diffs_mag[name][k][1])

                    before_val = float(fit_map[name][k][1])
                    denom = abs(before_val)
                    if denom < eps:
                        fatal("Cannot compute ratio for {} at index {}: abs(before value) is {:.17g} (too close to zero).".format(name, k, denom))
                    #endif

                    ratio = sysmag / denom
                    ratio_list.append([tmean, ratio])
                #endfor

                rhs_ratio = format_mathematica_list_pair(ratio_list)
                f.write(name + "_ratio = " + rhs_ratio + ";\n")
            #endfor
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
# Plotting (input series)
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
    for k in range(len(raw)):
        ent = raw[k]
        if (not isinstance(ent, (list, tuple))) or (len(ent) not in [2, 3]):
            fatal("Systematics series '{}' entry {} must be [tmean, sys] or [tmean, *, sys].".format(sys_name, k))
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

def save_input_asymmetry_canvases(fit_map, out_dir):
    np, plt = import_plot_deps()

    capsize = 3
    ms = 5.0
    xlim_t = (0.0, 1.30)
    x_label = r"$-t'\ (\mathrm{GeV}^{2})$"

    ylim_single = (-0.40, 0.40)
    ylim_double = (-0.80, 0.80)

    suffix_lu   = "GEchi2FitsALUsinphi"
    suffix_ul1  = "GEchi2FitsAULsinphi"
    suffix_ul2  = "GEchi2FitsAULsin2phi"
    suffix_ll0  = "GEchi2FitsALL"
    suffix_ll1  = "GEchi2FitsALLcosphi"

    groups = ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]

    for g in groups:
        v_lu  = g + suffix_lu
        v_ul1 = g + suffix_ul1
        v_ul2 = g + suffix_ul2
        v_ll0 = g + suffix_ll0
        v_ll1 = g + suffix_ll1

        for v in [v_lu, v_ul1, v_ul2, v_ll0, v_ll1]:
            if v not in fit_map:
                fatal("Cannot make plots: missing required input series '{}' in fit-results file.".format(v))
            #endif
        #endfor

        x_lu,  y_lu,  e_lu  = to_series(fit_map[v_lu],  np=np)
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

        sys_lu  = get_sys_band(fit_map, v_lu,  x_lu, np)
        sys_ul1 = get_sys_band(fit_map, v_ul1, x_lu, np)
        sys_ul2 = get_sys_band(fit_map, v_ul2, x_lu, np)
        sys_ll0 = get_sys_band(fit_map, v_ll0, x_lu, np)
        sys_ll1 = get_sys_band(fit_map, v_ll1, x_lu, np)

        sys_band_ul_panel = sys_ul1
        sys_band_ll_panel = sys_ll0

        if (sys_lu is None) and (sys_band_ul_panel is None) and (sys_band_ll_panel is None):
            warn("No systematics series found for bin '{}'; sys rectangles will be skipped.".format(g))
        #endif

        fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4))

        ax_left  = axes[0]
        ax_mid   = axes[1]
        ax_right = axes[2]

        draw_sys_bars(ax_left, x_lu, sys_lu, widths)
        draw_points(
            ax_left, x_lu, y_lu, e_lu,
            label=r"$F_{LU}^{\sin\phi}/F_{UU}$",
            color="tab:blue", marker="o", capsize=capsize, ms=ms
        )
        ax_left.set(xlim=xlim_t, ylim=ylim_single, xlabel=x_label, ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$")
        ax_left.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
        ax_left.grid(True, linestyle="--", alpha=0.6)

        draw_sys_bars(ax_mid, x_ul1, sys_band_ul_panel, widths)
        draw_points(
            ax_mid, x_ul1, y_ul1, e_ul1,
            label=r"$F_{UL}^{\sin\phi}/F_{UU}$",
            color="tab:red", marker="s", capsize=capsize, ms=ms
        )
        draw_points(
            ax_mid, x_ul2, y_ul2, e_ul2,
            label=r"$F_{UL}^{\sin2\phi}/F_{UU}$",
            color="tab:green", marker="^", capsize=capsize, ms=ms
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
            color="tab:purple", marker="o", capsize=capsize, ms=ms
        )
        draw_points(
            ax_right, x_ll1, y_ll1, e_ll1,
            label=r"$F_{LL}^{\cos\phi}/F_{UU}$",
            color="tab:orange", marker="s", capsize=capsize, ms=ms
        )
        ax_right.set(xlim=xlim_t, ylim=ylim_double, xlabel=x_label, ylabel=r"$F_{LL}^{(0,\cos\phi)}/F_{UU}$")
        ax_right.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
        ax_right.grid(True, linestyle="--", alpha=0.6)
        leg_right = ax_right.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_right.get_frame().set_alpha(0.9)

        plt.suptitle(r"$ep \rightarrow en\pi^{+}$" + " - " + xb_label(g), fontsize=14, y=0.98)
        plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])

        out_pdf = os.path.join(out_dir, "input_asymmetries_{}.pdf".format(g))
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

# =========================
# Plotting (migration sys vs -t')
# =========================

def save_migration_sys_canvases(fit_map, out_dir):
    np, plt = import_plot_deps()

    xlim_t = (0.0, 1.30)
    x_label = r"$-t'\ (\mathrm{GeV}^{2})$"

    suffix_ul1  = "GEchi2FitsAULsinphi"
    suffix_ul2  = "GEchi2FitsAULsin2phi"
    suffix_ll0  = "GEchi2FitsALL"
    suffix_ll1  = "GEchi2FitsALLcosphi"

    groups = ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]

    for g in groups:
        v_ul1 = g + suffix_ul1
        v_ul2 = g + suffix_ul2
        v_ll0 = g + suffix_ll0
        v_ll1 = g + suffix_ll1

        # These are injected by inject_migration_sys_for_plotting()
        v_ul1_sys = v_ul1 + "Sys"
        v_ul2_sys = v_ul2 + "Sys"
        v_ll0_sys = v_ll0 + "Sys"
        v_ll1_sys = v_ll1 + "Sys"

        for v in [v_ul1, v_ul2, v_ll0, v_ll1]:
            if v not in fit_map:
                fatal("Cannot make migration-sys plots: missing required input series '{}'.".format(v))
            #endif
        #endfor

        for v in [v_ul1_sys, v_ul2_sys, v_ll0_sys, v_ll1_sys]:
            if v not in fit_map:
                fatal("Cannot make migration-sys plots: missing injected sys series '{}'.".format(v))
            #endif
        #endfor

        # x reference from any base series (triples)
        x_ref, y_dummy, e_dummy = to_series(fit_map[v_ul1], np=np)

        sys_ul1 = get_sys_band(fit_map, v_ul1, x_ref, np)
        sys_ul2 = get_sys_band(fit_map, v_ul2, x_ref, np)
        sys_ll0 = get_sys_band(fit_map, v_ll0, x_ref, np)
        sys_ll1 = get_sys_band(fit_map, v_ll1, x_ref, np)

        if (sys_ul1 is None) or (sys_ul2 is None) or (sys_ll0 is None) or (sys_ll1 is None):
            fatal("Internal error: expected all sys series to exist for '{}' but at least one is missing.".format(g))
        #endif

        fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4))

        ax_ul = axes[0]
        ax_ll = axes[1]

        # UL subplot (two harmonics)
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

        # LL subplot (two harmonics)
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

        plt.suptitle(r"$ep \rightarrow en\pi^{+}$" + " - " + xb_label(g) + " - migration systematics", fontsize=14, y=0.98)
        plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])

        out_pdf = os.path.join(out_dir, "migration_systematics_{}.pdf".format(g))
        plt.savefig(out_pdf)
        plt.close(fig)

        sys.stdout.write("Saved plot: {}\n".format(out_pdf))
    #endfor
#enddef

def main():
    ap = argparse.ArgumentParser(
        description="Apply a 24x24 bin-migration matrix to selected enpi* GEchi2Fits* asymmetry lists."
    )
    ap.add_argument("migration_csv", help="CSV with 24 rows (bins 0-23) and 24 fraction columns per row.")
    ap.add_argument("fit_in_txt", help="Input Mathematica-style fit-results text file.")
    ap.add_argument("fit_out_name", help="Output file name only (no directory). Saved under output/rgc_enpi+_bin_migration_study/")

    ap.add_argument("--out-diff", default=None, help="Diff output file name only (no directory). Default: <fit_out_name with _diff inserted before extension>")
    ap.add_argument("--renormalize", action="store_true", help="Renormalize each migration row to sum to 1.0 before weighting.")
    ap.add_argument("--tol", type=float, default=1.0e-3, help="Tolerance for row-sum checks against 1.0 (strict mode).")

    args = ap.parse_args()

    out_dir = "output/rgc_enpi+_bin_migration_study"
    ensure_output_dir(out_dir)

    ensure_is_basename(args.fit_out_name, "fit_out_name")

    out_diff_name = args.out_diff
    if out_diff_name is None:
        out_diff_name = default_diff_basename(args.fit_out_name)
    #endif
    ensure_is_basename(out_diff_name, "out-diff")

    out_path = os.path.join(out_dir, args.fit_out_name)
    out_diff_path = os.path.join(out_dir, out_diff_name)

    weights, meta = read_migration_csv(args.migration_csv)
    fit_map = parse_fit_results_text(args.fit_in_txt)

    required = build_required_varnames()
    validate_fit_data(fit_map, required)

    migrated, diffs_mag, diffs_signed = compute_migrated_values(
        weights=weights,
        fit_map=fit_map,
        required_varnames=required,
        renormalize=args.renormalize,
        tol=args.tol
    )

    write_output_files(out_path, out_diff_path, migrated, diffs_mag, diffs_signed, fit_map)

    sys.stdout.write("Wrote migrated fit file: {}\n".format(out_path))
    sys.stdout.write("Wrote difference file:   {}\n".format(out_diff_path))

    fit_map_for_plots = inject_migration_sys_for_plotting(fit_map, diffs_mag, required)
    save_input_asymmetry_canvases(fit_map_for_plots, out_dir)

    # New: plots of migration systematic magnitudes vs -t' (both harmonics)
    save_migration_sys_canvases(fit_map_for_plots, out_dir)
#enddef

if __name__ == "__main__":
    main()
#endif