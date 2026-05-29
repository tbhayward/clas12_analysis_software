#!/usr/bin/env python3
# apply_bin_migration_to_asymmetries.py
#
# Usage:
# python3 apply_bin_migration_to_asymmetries.py mc_bin_migration.csv fit_results.txt fit_results_migrated.txt
#
# Purpose:
#   Apply the standard non-Fermi kinematic bin-migration matrix to the fitted
#   asymmetry points and treat the migration effect as a SIGNED correction.
#
# Inputs:
#   (1) Migration CSV with 24x24 per-bin composition/migration fractions.
#       Format:
#         #binNum,xBmin,xBmax,-tmin,-tmax,numEvents,fracfromMCbins0-23
#         0,0.1,0.25,0.05,0.25,nev,frac0,...,frac23
#         ...
#
#   (2) Fit-results text file in Mathematica-style assignments:
#         name = {{tmean, value, stat}, ...};
#       for 4 xB groups x 5 SF suffixes = 20 variables total, each of length 6.
#
# Outputs written under:
#   output/rgc_enpi+_bin_migration_study/
#
#   (1) Migrated/corrected fit file:
#         name = {{tmean_ALUsinphi, corrected_value, stat_original}, ...};
#
#   (2) Diff/correction/systematics file:
#         name_delta = {{tmean_ALUsinphi, signed_scaled_smoothed_delta}, ...};
#         name_abs   = {{tmean_ALUsinphi, scaled_smoothed_abs_delta}, ...};
#         name_ratio = {{tmean_ALUsinphi, signed_scaled_smoothed_delta/original}, ...};
#
#      Audit/debug quantities are also written:
#         name_raw_delta = {{tmean_ALUsinphi, raw_migrated-original}, ...};
#         name_raw_abs   = {{tmean_ALUsinphi, abs(raw_migrated-original)}, ...};
#         name_raw_ratio = {{tmean_ALUsinphi, raw_delta/original}, ...};
#
# IMPORTANT:
#   - The raw migration difference is:
#       raw_delta = raw_migrated - original
#
#   - The migration systematic magnitude is smoothed as a function of -tprime
#     within each xB bin and each modulation.
#
#   - The smoothed migration systematic is scaled by:
#       MIGRATION_SYS_SCALE = 1.5
#
#   - The exported signed delta uses:
#       final_delta = sign(raw_delta) * 1.5 * smooth(abs(raw_delta))
#
#   - Therefore the exported migrated/corrected fit file is:
#       corrected = original + final_delta
#
#   - The audit quantities preserve the raw matrix result before smoothing and
#     before the 1.5 scale factor.
#
# Important tmean convention:
#   For each xB group, this script uses the tmean values from:
#       <xB group>GEchi2FitsALUsinphi
#   as the canonical tmean values for ALL five fitted structure-function ratios.
#
#   This intentionally ignores small modulation-to-modulation tmean differences.
#   The migration calculation itself remains index-based.
#
# Plots:
#   (A) 4 PDF canvases for baseline/original points with migration sys rectangles:
#       polarized_structure_functions_baseline_<xBgroup>.pdf
#
#   (B) 4 PDF canvases for migrated/corrected points with migration sys rectangles:
#       polarized_structure_functions_migrated_<xBgroup>.pdf
#
#   (C) 4 PDF canvases showing only migration systematic magnitudes:
#       migration_systematics_<xBgroup>.pdf
#
# Plot requirements:
#   - Y ranges:
#       LU, UL: [-0.4, 0.4]
#       LL:     [-0.8, 0.8]
#   - No horizontal error bars.
#   - Systematic uncertainty drawn as symmetric gray rectangles from y=0 to +/-sys.
#   - UL panel uses the AUL sin(phi) migration systematic.
#   - LL panel uses the ALL migration systematic.

import argparse
import ast
import csv
import os
import re
import sys


MIGRATION_SYS_SCALE = 1.5

SMOOTHING_ENABLED = True

# Weak three-point smoothing kernel for the interior points:
#   smooth[i] = 0.25*y[i-1] + 0.50*y[i] + 0.25*y[i+1]
SMOOTHING_LEFT_WEIGHT = 0.25
SMOOTHING_CENTER_WEIGHT = 0.50
SMOOTHING_RIGHT_WEIGHT = 0.25

# One-sided edge smoothing:
#   smooth[0]  = 0.75*y[0]  + 0.25*y[1]
#   smooth[-1] = 0.75*y[-1] + 0.25*y[-2]
SMOOTHING_EDGE_SELF_WEIGHT = 0.75
SMOOTHING_EDGE_NEIGHBOR_WEIGHT = 0.25


# ---------------------------------------------------------------------------
# Basic helpers
# ---------------------------------------------------------------------------

def fatal(msg):
    sys.stderr.write("FATAL: " + msg + "\n")
    sys.exit(1)
#enddef


def warn(msg):
    sys.stderr.write("WARNING: " + msg + "\n")
#enddef


def info(msg):
    sys.stdout.write("INFO: " + msg + "\n")
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
        fatal("{} must be a file name only, no directory. Got: '{}'".format(what, path_like))
    #endif
#enddef


def sign_of(x):
    if x > 0.0:
        return 1.0
    #endif

    if x < 0.0:
        return -1.0
    #endif

    return 0.0
#enddef


# ---------------------------------------------------------------------------
# Naming conventions
# ---------------------------------------------------------------------------

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


def xb_group_names():
    return ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]
#enddef


def sf_suffixes():
    return [
        "GEchi2FitsALUsinphi",
        "GEchi2FitsAULsinphi",
        "GEchi2FitsAULsin2phi",
        "GEchi2FitsALL",
        "GEchi2FitsALLcosphi",
    ]
#enddef


def build_required_varnames():
    varnames = []

    for g in xb_group_names():
        for s in sf_suffixes():
            varnames.append(g + s)
        #endfor
    #endfor

    return varnames
#enddef


def global_bin_to_group_and_tindex(b):
    if b < 0 or b > 23:
        fatal("Global bin out of range: {}".format(b))
    #endif

    return (b // 6), (b % 6)
#enddef


def canonical_tmean_varname(xb_name):
    return xb_name + "GEchi2FitsALUsinphi"
#enddef


def get_canonical_tmeans_for_group(fit_map, xb_name):
    """
    Return the six canonical tmean values for this xB group.

    The canonical values are always taken from:
        <xb_name>GEchi2FitsALUsinphi

    These values are then used for every fitted structure-function ratio in this
    xB group in the migrated output, diff output, systematic-injected series,
    and diagnostic plots.
    """
    ref_name = canonical_tmean_varname(xb_name)

    if ref_name not in fit_map:
        fatal("Missing canonical tmean reference variable '{}'.".format(ref_name))
    #endif

    arr = fit_map[ref_name]
    if (not isinstance(arr, list)) or (len(arr) != 6):
        fatal("Canonical tmean reference '{}' must be a list of length 6.".format(ref_name))
    #endif

    tmeans = []
    for i in range(6):
        ent = arr[i]
        if (not isinstance(ent, (list, tuple))) or (len(ent) != 3):
            fatal("Canonical tmean reference '{}' entry {} must be [tmean, value, stat].".format(ref_name, i))
        #endif
        tmeans.append(float(ent[0]))
    #endfor

    return tmeans
#enddef


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def read_migration_csv(csv_path):
    """
    Migration CSV format per row:
      bin, xbmin, xbmax, tmin_n, tmax_n, nev, frac0..frac23

    Expected:
      24 rows
      30 columns per row = 6 metadata columns + 24 fractions

    Returns:
      weights: list[list[float]] shape (24, 24)
        row i = reconstructed target bin i
        column j = generated source bin j

      meta: list[(bin, xbmin, xbmax, tmin_n, tmax_n, nev)]
    """
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
                fatal(
                    "Migration CSV row {} has {} columns, expected 30 "
                    "(6 meta + 24 fractions).".format(row_idx, len(row))
                )
            #endif

            try:
                bin_num = int(float(row[0]))
                xbmin = float(row[1])
                xbmax = float(row[2])
                tmin_n = float(row[3])
                tmax_n = float(row[4])
                nev = int(float(row[5]))
            except Exception as e:
                fatal("Failed parsing meta columns on CSV row {}: {}".format(row_idx, str(e)))
            #endtry

            fracs = []
            for j in range(24):
                try:
                    fracs.append(float(row[6 + j]))
                except Exception as e:
                    fatal("Failed parsing fraction column j={} on CSV row {}: {}".format(j, row_idx, str(e)))
                #endtry
            #endfor

            weights.append(fracs)
            meta.append((bin_num, xbmin, xbmax, tmin_n, tmax_n, nev))
        #endfor
    #endwith

    if len(weights) != 24:
        fatal("Migration CSV has {} data rows; expected 24 bins.".format(len(weights)))
    #endif

    for i in range(24):
        if meta[i][0] != i:
            fatal("Migration CSV binNum mismatch: row {} has binNum {}, expected {}.".format(i, meta[i][0], i))
        #endif
    #endfor

    return weights, meta
#enddef


def parse_fit_results_text(txt_path):
    if not os.path.isfile(txt_path):
        fatal("Fit-results text file not found: " + txt_path)
    #endif

    with open(txt_path, "r") as f:
        text = f.read()
    #endwith

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
            fatal("Failed to parse assignment for variable '{}' in '{}': {}".format(name, txt_path, str(e)))
        #endtry

        out[name] = val
    #endfor

    return out
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
            fatal(
                "Variable '{}' must be a list of length 6; got type {} length {}.".format(
                    v,
                    type(arr).__name__,
                    len(arr) if isinstance(arr, list) else -1,
                )
            )
        #endif

        for k in range(6):
            ent = arr[k]
            if (not isinstance(ent, (list, tuple))) or (len(ent) != 3):
                fatal("Variable '{}' entry {} must be [tmean, value, stat].".format(v, k))
            #endif

            try:
                float(ent[0])
                float(ent[1])
                float(ent[2])
            except Exception as e:
                fatal("Variable '{}' entry {} contains non-numeric value: {}".format(v, k, str(e)))
            #endtry
        #endfor
    #endfor

    for g in xb_group_names():
        get_canonical_tmeans_for_group(fit_map, g)
    #endfor
#enddef


# ---------------------------------------------------------------------------
# Smoothing and scale-factor treatment
# ---------------------------------------------------------------------------

def smooth_abs_values(abs_values):
    """
    Apply a weak three-point smoothing filter to a one-dimensional six-point
    migration-systematic magnitude sequence.

    Interior points:
        y_smooth[i] = 0.25*y[i-1] + 0.50*y[i] + 0.25*y[i+1]

    Edge points:
        y_smooth[0]  = 0.75*y[0]  + 0.25*y[1]
        y_smooth[-1] = 0.75*y[-1] + 0.25*y[-2]

    The function is intentionally weak: it reduces isolated bin-to-bin
    fluctuations while preserving broad trends in -tprime.
    """
    n = len(abs_values)

    if n == 0:
        return []
    #endif

    if n == 1:
        return [float(abs_values[0])]
    #endif

    smoothed = [0.0 for _ in range(n)]

    smoothed[0] = (
        SMOOTHING_EDGE_SELF_WEIGHT * float(abs_values[0])
        + SMOOTHING_EDGE_NEIGHBOR_WEIGHT * float(abs_values[1])
    )

    for i in range(1, n - 1):
        smoothed[i] = (
            SMOOTHING_LEFT_WEIGHT * float(abs_values[i - 1])
            + SMOOTHING_CENTER_WEIGHT * float(abs_values[i])
            + SMOOTHING_RIGHT_WEIGHT * float(abs_values[i + 1])
        )
    #endfor

    smoothed[n - 1] = (
        SMOOTHING_EDGE_SELF_WEIGHT * float(abs_values[n - 1])
        + SMOOTHING_EDGE_NEIGHBOR_WEIGHT * float(abs_values[n - 2])
    )

    return smoothed
#enddef


def apply_smoothing_and_scale_to_diffs(raw_diffs_signed, raw_diffs_mag):
    """
    Convert raw migration deltas to final exported deltas and systematics.

    Input:
      raw_diffs_signed[var] = [[tmean, raw_delta], ...]
      raw_diffs_mag[var]    = [[tmean, abs(raw_delta)], ...]

    Output:
      final_diffs_signed[var] = [[tmean, sign(raw_delta)*1.5*smooth(abs(raw_delta))], ...]
      final_diffs_mag[var]    = [[tmean, 1.5*smooth(abs(raw_delta))], ...]

    This is done independently for every modulation and every xB group, because
    each variable has its own six-point -tprime dependence.
    """
    final_diffs_signed = {}
    final_diffs_mag = {}

    for name in raw_diffs_signed:
        if name not in raw_diffs_mag:
            fatal("Internal error: '{}' missing from raw_diffs_mag.".format(name))
        #endif

        signed_series = raw_diffs_signed[name]
        mag_series = raw_diffs_mag[name]

        if len(signed_series) != len(mag_series):
            fatal("Internal error: length mismatch in raw diff series for '{}'.".format(name))
        #endif

        tmeans = []
        signs = []
        mags = []

        for i in range(len(signed_series)):
            t_signed = float(signed_series[i][0])
            d_raw = float(signed_series[i][1])
            t_mag = float(mag_series[i][0])
            m_raw = float(mag_series[i][1])

            if abs(t_signed - t_mag) > 1.0e-9:
                fatal("Internal error: tmean mismatch between raw signed and raw mag series for '{}' index {}.".format(name, i))
            #endif

            tmeans.append(t_signed)
            signs.append(sign_of(d_raw))
            mags.append(abs(m_raw))
        #endfor

        if SMOOTHING_ENABLED:
            mags_smoothed = smooth_abs_values(mags)
        else:
            mags_smoothed = list(mags)
        #endif

        final_signed_list = []
        final_mag_list = []

        for i in range(len(mags_smoothed)):
            final_mag = MIGRATION_SYS_SCALE * float(mags_smoothed[i])
            final_signed = signs[i] * final_mag

            final_signed_list.append([tmeans[i], final_signed])
            final_mag_list.append([tmeans[i], final_mag])
        #endfor

        final_diffs_signed[name] = final_signed_list
        final_diffs_mag[name] = final_mag_list
    #endfor

    return final_diffs_mag, final_diffs_signed
#enddef


def print_smoothing_and_scale_summary(raw_diffs_mag, final_diffs_mag):
    info("Migration systematic treatment:")
    info("  Raw migration differences are first converted to magnitudes |delta_raw|.")
    info("  A weak 3-point smoothing filter is applied independently to each modulation and xB bin.")
    info("  Interior smoothing: 0.25*left + 0.50*center + 0.25*right.")
    info("  Edge smoothing: 0.75*edge + 0.25*nearest_neighbor.")
    info("  The smoothed magnitudes are then multiplied by {:.3f}.".format(MIGRATION_SYS_SCALE))
    info("  The exported signed delta keeps the original raw sign but uses the scaled/smoothed magnitude.")
    info("  Raw audit series are also written as *_raw_delta, *_raw_abs, and *_raw_ratio.")

    max_raw = 0.0
    max_final = 0.0

    for name in raw_diffs_mag:
        for ent in raw_diffs_mag[name]:
            max_raw = max(max_raw, abs(float(ent[1])))
        #endfor
    #endfor

    for name in final_diffs_mag:
        for ent in final_diffs_mag[name]:
            max_final = max(max_final, abs(float(ent[1])))
        #endfor
    #endfor

    info("  Maximum raw migration magnitude over all points: {:.9f}".format(max_raw))
    info("  Maximum final scaled/smoothed magnitude over all points: {:.9f}".format(max_final))
#enddef


# ---------------------------------------------------------------------------
# Migration calculation
# ---------------------------------------------------------------------------

def compute_raw_migrated_values(weights, fit_map, renormalize, tol):
    """
    For each target reconstructed bin b:
      raw_migrated(b) = sum_j w_bj * value(j)

    where:
      b = target reconstructed bin
      j = generated/source bin

    Raw correction convention:
      raw_delta = raw_migrated - original

    Important:
      The tmean written to every output point is the canonical ALUsinphi tmean
      for that xB group and t index.
    """
    raw_migrated = {}
    raw_diffs_mag = {}
    raw_diffs_signed = {}

    for suffix in sf_suffixes():
        for xb_group_idx in range(4):
            xb_name = xb_group_name_from_index(xb_group_idx)
            varname = xb_name + suffix
            canonical_tmeans = get_canonical_tmeans_for_group(fit_map, xb_name)

            out_list = []
            diff_mag_list = []
            diff_signed_list = []

            for t_idx in range(6):
                target_global_bin = xb_group_idx * 6 + t_idx
                row = weights[target_global_bin]

                row_sum = 0.0
                for w in row:
                    row_sum += w
                #endfor

                if row_sum == 0.0:
                    fatal("Migration row sum is zero for target bin {}.".format(target_global_bin))
                #endif

                if (abs(row_sum - 1.0) > tol) and (not renormalize):
                    fatal(
                        "Migration row sum for target bin {} is {:.9f}, not within tol {:.9f} of 1.0.\n"
                        "If this is only rounding or an intentional choice, rerun with --renormalize.".format(
                            target_global_bin,
                            row_sum,
                            tol,
                        )
                    )
                #endif

                tmean_canonical = float(canonical_tmeans[t_idx])
                aval_target = float(fit_map[varname][t_idx][1])
                stat_target = float(fit_map[varname][t_idx][2])

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

                out_list.append([tmean_canonical, migrated_val, stat_target])
                diff_mag_list.append([tmean_canonical, delta_mag])
                diff_signed_list.append([tmean_canonical, delta_signed])
            #endfor

            raw_migrated[varname] = out_list
            raw_diffs_mag[varname] = diff_mag_list
            raw_diffs_signed[varname] = diff_signed_list
        #endfor
    #endfor

    return raw_migrated, raw_diffs_mag, raw_diffs_signed
#enddef


def build_final_corrected_values(fit_map, final_diffs_signed):
    """
    Build corrected output values using:
        corrected = original + final_delta

    where final_delta is the signed scaled/smoothed migration effect.

    Statistical uncertainties are kept unchanged from the original fit input.
    """
    corrected = {}

    for g in xb_group_names():
        canonical_tmeans = get_canonical_tmeans_for_group(fit_map, g)

        for s in sf_suffixes():
            name = g + s

            if name not in final_diffs_signed:
                fatal("Internal error: final_diffs_signed missing '{}'.".format(name))
            #endif

            out_list = []

            for i in range(6):
                tmean = float(canonical_tmeans[i])
                original_val = float(fit_map[name][i][1])
                stat = float(fit_map[name][i][2])
                delta = float(final_diffs_signed[name][i][1])

                corrected_val = original_val + delta
                out_list.append([tmean, corrected_val, stat])
            #endfor

            corrected[name] = out_list
        #endfor
    #endfor

    return corrected
#enddef


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

def format_mathematica_list_triple(lst):
    parts = []

    for trip in lst:
        parts.append(
            "{"
            + "{:.9f}".format(float(trip[0]))
            + ", "
            + "{:.9f}".format(float(trip[1]))
            + ", "
            + "{:.9f}".format(float(trip[2]))
            + "}"
        )
    #endfor

    return "{" + ", ".join(parts) + "}"
#enddef


def format_mathematica_list_pair(lst):
    parts = []

    for pair in lst:
        parts.append(
            "{"
            + "{:.9f}".format(float(pair[0]))
            + ", "
            + "{:.9f}".format(float(pair[1]))
            + "}"
        )
    #endfor

    return "{" + ", ".join(parts) + "}"
#enddef


def build_ratio_list(diff_signed, fit_map, name):
    eps = 1.0e-15
    ratio_list = []

    for k in range(6):
        tmean = float(diff_signed[name][k][0])
        delta = float(diff_signed[name][k][1])
        before_val = float(fit_map[name][k][1])

        if abs(before_val) < eps:
            fatal(
                "Cannot compute signed ratio for {} at index {}: "
                "baseline value is {:.17g}.".format(name, k, before_val)
            )
        #endif

        ratio = delta / before_val
        ratio_list.append([tmean, ratio])
    #endfor

    return ratio_list
#enddef


def write_output_files(
    out_path,
    out_diff_path,
    corrected,
    final_diffs_mag,
    final_diffs_signed,
    raw_diffs_mag,
    raw_diffs_signed,
    fit_map,
):
    """
    Writes:
      - out_path:
          corrected triples using final signed scaled/smoothed delta

      - out_diff_path:
          final signed delta, final absolute systematic, final signed ratio

        and audit quantities:
          raw signed delta, raw absolute delta, raw signed ratio
    """
    with open(out_path, "w") as f:
        for g in xb_group_names():
            for s in sf_suffixes():
                name = g + s

                if name not in corrected:
                    fatal("Internal error: corrected missing variable '{}'.".format(name))
                #endif

                rhs = format_mathematica_list_triple(corrected[name])
                f.write(name + " = " + rhs + ";\n")
            #endfor

            f.write("\n")
        #endfor
    #endwith

    with open(out_diff_path, "w") as f:
        f.write("# Migration systematic treatment applied in this file:\n")
        f.write("#   raw_delta = raw_migrated - original\n")
        f.write("#   abs(raw_delta) is weakly smoothed versus -tprime within each xB bin and modulation\n")
        f.write("#   final_abs = {:.9f} * smooth(abs(raw_delta))\n".format(MIGRATION_SYS_SCALE))
        f.write("#   final_delta = sign(raw_delta) * final_abs\n")
        f.write("#   *_delta, *_abs, *_ratio use the final scaled/smoothed values\n")
        f.write("#   *_raw_delta, *_raw_abs, *_raw_ratio preserve unsmoothed/unscaled audit values\n\n")

        for g in xb_group_names():
            for s in sf_suffixes():
                name = g + s

                if name not in final_diffs_signed:
                    fatal("Internal error: final_diffs_signed missing variable '{}'.".format(name))
                #endif

                if name not in final_diffs_mag:
                    fatal("Internal error: final_diffs_mag missing variable '{}'.".format(name))
                #endif

                if name not in raw_diffs_signed:
                    fatal("Internal error: raw_diffs_signed missing variable '{}'.".format(name))
                #endif

                if name not in raw_diffs_mag:
                    fatal("Internal error: raw_diffs_mag missing variable '{}'.".format(name))
                #endif

                if name not in fit_map:
                    fatal("Internal error: fit_map missing variable '{}' for ratio.".format(name))
                #endif

                rhs_delta = format_mathematica_list_pair(final_diffs_signed[name])
                f.write(name + "_delta = " + rhs_delta + ";\n")

                rhs_abs = format_mathematica_list_pair(final_diffs_mag[name])
                f.write(name + "_abs = " + rhs_abs + ";\n")

                ratio_list = build_ratio_list(final_diffs_signed, fit_map, name)
                rhs_ratio = format_mathematica_list_pair(ratio_list)
                f.write(name + "_ratio = " + rhs_ratio + ";\n")

                rhs_raw_delta = format_mathematica_list_pair(raw_diffs_signed[name])
                f.write(name + "_raw_delta = " + rhs_raw_delta + ";\n")

                rhs_raw_abs = format_mathematica_list_pair(raw_diffs_mag[name])
                f.write(name + "_raw_abs = " + rhs_raw_abs + ";\n")

                raw_ratio_list = build_ratio_list(raw_diffs_signed, fit_map, name)
                rhs_raw_ratio = format_mathematica_list_pair(raw_ratio_list)
                f.write(name + "_raw_ratio = " + rhs_raw_ratio + ";\n\n")
            #endfor
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


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def import_plot_deps():
    try:
        import numpy as np
    except Exception as e:
        fatal("numpy import failed, needed for plots: {}".format(str(e)))
    #endtry

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        fatal("matplotlib import failed, needed for plots: {}".format(str(e)))
    #endtry

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
        x = x[idx]
        y = y[idx]
        e = e[idx]
    #endif

    return x, y, e
#enddef


def to_sys_series(pairs, negate_x=True, sort=True, np=None):
    if np is None:
        fatal("Internal error: numpy not provided to to_sys_series().")
    #endif

    arr = np.array(pairs, dtype=float)

    x = -arr[:, 0] if negate_x else arr[:, 0]
    s = arr[:, 1]

    if sort:
        idx = np.argsort(x)
        x = x[idx]
        s = s[idx]
    #endif

    return x, s
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
        fatal(
            "Systematics series '{}' exists but has length {}, expected {}.".format(
                sys_name,
                len(raw) if isinstance(raw, list) else -1,
                len(x_ref),
            )
        )
    #endif

    pairs = []

    for k in range(len(raw)):
        ent = raw[k]

        if (not isinstance(ent, (list, tuple))) or (len(ent) not in [2, 3]):
            fatal("Systematics series '{}' entry {} must be [tmean, sys] or [tmean, dummy, sys].".format(sys_name, k))
        #endif

        tmean = float(ent[0])

        if len(ent) == 2:
            sysv = float(ent[1])
        else:
            sysv = float(ent[2])
        #endif

        pairs.append([tmean, abs(sysv)])
    #endfor

    x, s = to_sys_series(pairs, np=np)

    if (len(x) == len(x_ref)) and (not np.allclose(x, x_ref, rtol=0.0, atol=1.0e-6)):
        warn("x mismatch between '{}' systematics and stat series; using index pairing.".format(sys_name))
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

    ax.bar(
        x_centers,
        sys_band,
        width=widths,
        bottom=0.0,
        color="0.7",
        alpha=0.25,
        linewidth=0.0,
        align="center",
        zorder=0,
    )

    ax.bar(
        x_centers,
        -sys_band,
        width=widths,
        bottom=0.0,
        color="0.7",
        alpha=0.25,
        linewidth=0.0,
        align="center",
        zorder=0,
    )
#enddef


def draw_points(ax, x, y, e, label, color, marker, capsize, ms):
    ax.errorbar(
        x,
        y,
        yerr=e,
        fmt=marker,
        color=color,
        ecolor=color,
        capsize=capsize,
        markersize=ms,
        linestyle="None",
        label=label,
        zorder=2,
    )
#enddef


def xb_label(group_name):
    labels = {
        "enpiLowxB": r"$0.10 < x_{B} < 0.25$",
        "enpiMidLowxB": r"$0.25 < x_{B} < 0.35$",
        "enpiMidHighxB": r"$0.35 < x_{B} < 0.45$",
        "enpiHighxB": r"$0.45 < x_{B} < 0.60$",
    }

    return labels.get(group_name, group_name)
#enddef


# ---------------------------------------------------------------------------
# Plot data preparation
# ---------------------------------------------------------------------------

def make_canonical_plot_series_for_group(fit_map, g, suffix, x_ref, np):
    """
    Convert one fitted series to x,y,e for plotting, but force x to x_ref.

    The y and e arrays are still sorted according to the original x ordering,
    which should preserve the t-bin order. In practice the fit files are already
    in t-bin order, but this keeps behavior consistent with previous plotting.
    """
    varname = g + suffix
    x_raw, y, e = to_series(fit_map[varname], np=np)

    if len(x_raw) != len(x_ref):
        fatal("Series '{}' has {} points, expected {}.".format(varname, len(x_raw), len(x_ref)))
    #endif

    return x_ref.copy(), y, e
#enddef


def canonicalize_fit_map_tmeans(fit_map, required_varnames):
    """
    Return a copy of fit_map where every required variable uses the ALUsinphi
    tmean values for its xB group.

    This is useful for plotting and for injected systematics. It also keeps
    non-required variables unchanged.
    """
    out = dict(fit_map)

    for g in xb_group_names():
        canonical_tmeans = get_canonical_tmeans_for_group(fit_map, g)

        for s in sf_suffixes():
            name = g + s
            if name not in fit_map:
                continue
            #endif

            new_list = []

            for i in range(6):
                val = float(fit_map[name][i][1])
                stat = float(fit_map[name][i][2])
                new_list.append([float(canonical_tmeans[i]), val, stat])
            #endfor

            out[name] = new_list
        #endfor
    #endfor

    return out
#enddef


def inject_migration_sys_for_plotting(fit_map, final_diffs_mag, required_varnames):
    """
    Inject final scaled/smoothed migration systematic magnitudes as:
      <varname>Sys = {{canonical_tmean, final_abs_delta}, ...}
    """
    out = dict(fit_map)

    for v in required_varnames:
        if v not in final_diffs_mag:
            fatal("Internal error: final_diffs_mag missing '{}' for sys injection.".format(v))
        #endif

        sys_pairs = []

        for ent in final_diffs_mag[v]:
            if (not isinstance(ent, (list, tuple))) or (len(ent) != 2):
                fatal("Internal error: final_diffs_mag['{}'] entry must be [tmean, sysmag].".format(v))
            #endif

            tmean = float(ent[0])
            sysmag = abs(float(ent[1]))
            sys_pairs.append([tmean, sysmag])
        #endfor

        out[v + "Sys"] = sys_pairs
    #endfor

    return out
#enddef


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def save_polarized_structure_function_canvases(fit_map, out_dir, file_prefix):
    """
    Makes 1x3 canvases:
      - LU
      - UL: sin(phi), sin(2phi)
      - LL: constant, cos(phi)

    All x-values inside each xB group are forced to use the ALUsinphi x-axis.
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

    for g in xb_group_names():
        v_lu = g + suffix_lu
        v_ul1 = g + suffix_ul1
        v_ul2 = g + suffix_ul2
        v_ll0 = g + suffix_ll0
        v_ll1 = g + suffix_ll1

        for v in [v_lu, v_ul1, v_ul2, v_ll0, v_ll1]:
            if v not in fit_map:
                fatal("Cannot make plots: missing required input series '{}'.".format(v))
            #endif
        #endfor

        x_ref, y_lu, e_lu = to_series(fit_map[v_lu], np=np)

        x_lu = x_ref.copy()
        x_ul1, y_ul1, e_ul1 = make_canonical_plot_series_for_group(fit_map, g, suffix_ul1, x_ref, np)
        x_ul2, y_ul2, e_ul2 = make_canonical_plot_series_for_group(fit_map, g, suffix_ul2, x_ref, np)
        x_ll0, y_ll0, e_ll0 = make_canonical_plot_series_for_group(fit_map, g, suffix_ll0, x_ref, np)
        x_ll1, y_ll1, e_ll1 = make_canonical_plot_series_for_group(fit_map, g, suffix_ll1, x_ref, np)

        left, right, widths = compute_bin_edges_from_centers(x_ref, np)

        sys_lu = get_sys_band(fit_map, v_lu, x_ref, np)
        sys_ul1 = get_sys_band(fit_map, v_ul1, x_ref, np)
        sys_ul2 = get_sys_band(fit_map, v_ul2, x_ref, np)
        sys_ll0 = get_sys_band(fit_map, v_ll0, x_ref, np)
        sys_ll1 = get_sys_band(fit_map, v_ll1, x_ref, np)

        sys_band_ul_panel = sys_ul1
        sys_band_ll_panel = sys_ll0

        fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4))

        ax_left = axes[0]
        ax_mid = axes[1]
        ax_right = axes[2]

        draw_sys_bars(ax_left, x_lu, sys_lu, widths)
        draw_points(
            ax_left,
            x_lu,
            y_lu,
            e_lu,
            label=r"$F_{LU}^{\sin\phi}/F_{UU}$",
            color="tab:blue",
            marker="o",
            capsize=3,
            ms=5.0,
        )
        ax_left.set(
            xlim=xlim_t,
            ylim=ylim_single,
            xlabel=x_label,
            ylabel=r"$F_{LU}^{\sin\phi}/F_{UU}$",
        )
        ax_left.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
        ax_left.grid(True, linestyle="--", alpha=0.6)

        draw_sys_bars(ax_mid, x_ref, sys_band_ul_panel, widths)
        draw_points(
            ax_mid,
            x_ul1,
            y_ul1,
            e_ul1,
            label=r"$F_{UL}^{\sin\phi}/F_{UU}$",
            color="tab:red",
            marker="s",
            capsize=3,
            ms=5.0,
        )
        draw_points(
            ax_mid,
            x_ul2,
            y_ul2,
            e_ul2,
            label=r"$F_{UL}^{\sin2\phi}/F_{UU}$",
            color="tab:green",
            marker="^",
            capsize=3,
            ms=5.0,
        )
        ax_mid.set(
            xlim=xlim_t,
            ylim=ylim_single,
            xlabel=x_label,
            ylabel=r"$F_{UL}^{\sin(n\phi)}/F_{UU}$",
        )
        ax_mid.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
        ax_mid.grid(True, linestyle="--", alpha=0.6)
        leg_mid = ax_mid.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_mid.get_frame().set_alpha(0.9)

        draw_sys_bars(ax_right, x_ref, sys_band_ll_panel, widths)
        draw_points(
            ax_right,
            x_ll0,
            y_ll0,
            e_ll0,
            label=r"$F_{LL}/F_{UU}$",
            color="tab:purple",
            marker="o",
            capsize=3,
            ms=5.0,
        )
        draw_points(
            ax_right,
            x_ll1,
            y_ll1,
            e_ll1,
            label=r"$F_{LL}^{\cos\phi}/F_{UU}$",
            color="tab:orange",
            marker="s",
            capsize=3,
            ms=5.0,
        )
        ax_right.set(
            xlim=xlim_t,
            ylim=ylim_double,
            xlabel=x_label,
            ylabel=r"$F_{LL}^{(0,\cos\phi)}/F_{UU}$",
        )
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


def save_migration_sys_canvases(fit_map, out_dir):
    """
    Make 4 PDFs, one per xB group, showing only final scaled/smoothed migration
    systematic magnitudes.

    All x-values are taken from the ALUsinphi canonical x-axis in each group.
    """
    np, plt = import_plot_deps()

    xlim_t = (0.0, 1.30)
    x_label = r"$-t'\ (\mathrm{GeV}^{2})$"

    suffix_lu = "GEchi2FitsALUsinphi"
    suffix_ul1 = "GEchi2FitsAULsinphi"
    suffix_ul2 = "GEchi2FitsAULsin2phi"
    suffix_ll0 = "GEchi2FitsALL"
    suffix_ll1 = "GEchi2FitsALLcosphi"

    for g in xb_group_names():
        v_lu = g + suffix_lu
        v_ul1 = g + suffix_ul1
        v_ul2 = g + suffix_ul2
        v_ll0 = g + suffix_ll0
        v_ll1 = g + suffix_ll1

        for v in [v_lu, v_ul1, v_ul2, v_ll0, v_ll1]:
            if v not in fit_map:
                fatal("Cannot make migration-sys plots: missing required input series '{}'.".format(v))
            #endif
        #endfor

        x_ref, y_dummy, e_dummy = to_series(fit_map[v_lu], np=np)

        sys_lu = get_sys_band(fit_map, v_lu, x_ref, np)
        sys_ul1 = get_sys_band(fit_map, v_ul1, x_ref, np)
        sys_ul2 = get_sys_band(fit_map, v_ul2, x_ref, np)
        sys_ll0 = get_sys_band(fit_map, v_ll0, x_ref, np)
        sys_ll1 = get_sys_band(fit_map, v_ll1, x_ref, np)

        if (sys_lu is None) or (sys_ul1 is None) or (sys_ul2 is None) or (sys_ll0 is None) or (sys_ll1 is None):
            fatal("Internal error: expected all sys series for '{}' but at least one is missing.".format(g))
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
        ax_lu.set(
            xlim=xlim_t,
            ylim=(0.0, 1.2 * y_lu_max),
            xlabel=x_label,
            ylabel="Migration systematic magnitude",
        )
        ax_lu.grid(True, linestyle="--", alpha=0.6)
        leg_lu = ax_lu.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_lu.get_frame().set_alpha(0.9)

        y_ul_max = float(np.max(np.array([np.max(sys_ul1), np.max(sys_ul2)], dtype=float)))
        if y_ul_max <= 0.0:
            y_ul_max = 1.0e-6
        #endif

        ax_ul.plot(x_ref, sys_ul1, marker="o", linestyle="None", label=r"$|\Delta(F_{UL}^{\sin\phi}/F_{UU})|$")
        ax_ul.plot(x_ref, sys_ul2, marker="s", linestyle="None", label=r"$|\Delta(F_{UL}^{\sin2\phi}/F_{UU})|$")
        ax_ul.set(
            xlim=xlim_t,
            ylim=(0.0, 1.2 * y_ul_max),
            xlabel=x_label,
            ylabel="Migration systematic magnitude",
        )
        ax_ul.grid(True, linestyle="--", alpha=0.6)
        leg_ul = ax_ul.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_ul.get_frame().set_alpha(0.9)

        y_ll_max = float(np.max(np.array([np.max(sys_ll0), np.max(sys_ll1)], dtype=float)))
        if y_ll_max <= 0.0:
            y_ll_max = 1.0e-6
        #endif

        ax_ll.plot(x_ref, sys_ll0, marker="o", linestyle="None", label=r"$|\Delta(F_{LL}/F_{UU})|$")
        ax_ll.plot(x_ref, sys_ll1, marker="s", linestyle="None", label=r"$|\Delta(F_{LL}^{\cos\phi}/F_{UU})|$")
        ax_ll.set(
            xlim=xlim_t,
            ylim=(0.0, 1.2 * y_ll_max),
            xlabel=x_label,
            ylabel="Migration systematic magnitude",
        )
        ax_ll.grid(True, linestyle="--", alpha=0.6)
        leg_ll = ax_ll.legend(frameon=True, edgecolor="black", fontsize=10, loc="upper right")
        leg_ll.get_frame().set_alpha(0.9)

        plt.suptitle(r"$ep \rightarrow en\pi^{+}$" + " - " + xb_label(g) + " - bin-migration systematics", fontsize=14, y=0.98)
        plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])

        out_pdf = os.path.join(out_dir, "migration_systematics_{}.pdf".format(g))
        plt.savefig(out_pdf)
        plt.close(fig)

        sys.stdout.write("Saved plot: {}\n".format(out_pdf))
    #endfor
#enddef


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Apply standard 24x24 bin-migration mixture to enpi* GEchi2Fits* asymmetry lists."
    )

    ap.add_argument(
        "migration_csv",
        help="Bin migration CSV with 24 rows and 30 columns: 6 metadata + 24 fractions.",
    )

    ap.add_argument(
        "fit_txt",
        help="Input Mathematica-style fit-results text file.",
    )

    ap.add_argument(
        "fit_out_name",
        help="Output file name only, no directory. Written under output/rgc_enpi+_bin_migration_study/",
    )

    ap.add_argument(
        "--out-diff",
        default=None,
        help="Diff output file name only. Default: <fit_out_name with _diff inserted>.",
    )

    ap.add_argument(
        "--renormalize",
        action="store_true",
        help="Renormalize each migration row to sum to 1.0 after reading.",
    )

    ap.add_argument(
        "--tol",
        type=float,
        default=1.0e-3,
        help="Tolerance for row-sum checks against 1.0 in strict mode.",
    )

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

    info("Starting bin-migration systematic calculation.")
    info("Input migration CSV: {}".format(args.migration_csv))
    info("Input fit-results file: {}".format(args.fit_txt))
    info("Output migrated/corrected file: {}".format(out_path))
    info("Output diff/systematics file: {}".format(out_diff_path))

    warn("Reviewer-coordinated migration systematic scale factor is active: {:.3f}".format(MIGRATION_SYS_SCALE))

    if SMOOTHING_ENABLED:
        warn("Weak smoothing of migration systematic magnitudes is active.")
        warn("Smoothing is applied to |raw_delta| versus -tprime for each modulation in each xB bin.")
        warn("Interior kernel: 0.25*left + 0.50*center + 0.25*right.")
        warn("Edge kernel: 0.75*edge + 0.25*nearest_neighbor.")
    else:
        warn("Smoothing is disabled.")
    #endif

    weights, meta = read_migration_csv(args.migration_csv)

    fit_map_original = parse_fit_results_text(args.fit_txt)
    required = build_required_varnames()
    validate_fit_data(fit_map_original, required)

    fit_map_canonical = canonicalize_fit_map_tmeans(fit_map_original, required)

    raw_migrated, raw_diffs_mag, raw_diffs_signed = compute_raw_migrated_values(
        weights=weights,
        fit_map=fit_map_canonical,
        renormalize=args.renormalize,
        tol=args.tol,
    )

    final_diffs_mag, final_diffs_signed = apply_smoothing_and_scale_to_diffs(
        raw_diffs_signed=raw_diffs_signed,
        raw_diffs_mag=raw_diffs_mag,
    )

    print_smoothing_and_scale_summary(
        raw_diffs_mag=raw_diffs_mag,
        final_diffs_mag=final_diffs_mag,
    )

    corrected = build_final_corrected_values(
        fit_map=fit_map_canonical,
        final_diffs_signed=final_diffs_signed,
    )

    write_output_files(
        out_path=out_path,
        out_diff_path=out_diff_path,
        corrected=corrected,
        final_diffs_mag=final_diffs_mag,
        final_diffs_signed=final_diffs_signed,
        raw_diffs_mag=raw_diffs_mag,
        raw_diffs_signed=raw_diffs_signed,
        fit_map=fit_map_canonical,
    )

    sys.stdout.write("Wrote migrated/corrected fit file: {}\n".format(out_path))
    sys.stdout.write("Wrote difference/systematics file: {}\n".format(out_diff_path))

    baseline_for_plots = inject_migration_sys_for_plotting(fit_map_canonical, final_diffs_mag, required)
    corrected_for_plots = inject_migration_sys_for_plotting(corrected, final_diffs_mag, required)

    save_polarized_structure_function_canvases(
        baseline_for_plots,
        out_dir,
        file_prefix="polarized_structure_functions_baseline",
    )

    save_polarized_structure_function_canvases(
        corrected_for_plots,
        out_dir,
        file_prefix="polarized_structure_functions_migrated",
    )

    save_migration_sys_canvases(baseline_for_plots, out_dir)

    sys.stdout.write("\nSign convention reminder:\n")
    sys.stdout.write("  raw_delta = raw_migrated - original\n")
    sys.stdout.write("  final_abs = {:.3f} * smooth(abs(raw_delta))\n".format(MIGRATION_SYS_SCALE))
    sys.stdout.write("  final_delta = sign(raw_delta) * final_abs\n")
    sys.stdout.write("  corrected value written to migrated/corrected file is original + final_delta.\n")
    sys.stdout.write("  final_delta > 0 shifts asymmetry upward, numerically larger.\n")

    sys.stdout.write("\nAudit reminder:\n")
    sys.stdout.write("  *_delta, *_abs, *_ratio use scaled/smoothed migration systematics.\n")
    sys.stdout.write("  *_raw_delta, *_raw_abs, *_raw_ratio preserve the unsmoothed/unscaled matrix result.\n")

    sys.stdout.write("\ntmean convention reminder:\n")
    sys.stdout.write("  All output tmean values use the ALUsinphi tmean values for the corresponding xB group.\n")
#enddef


if __name__ == "__main__":
    main()
#endif