#!/usr/bin/env python3
# apply_bin_migration_to_asymmetries.py
#
# Usage (tcsh one-liner example):
# python3 apply_bin_migration_to_asymmetries.py mc_bin_migration.csv fit_results.txt fit_results_migrated.txt

import argparse
import ast
import csv
import os
import re
import sys

def fatal(msg):
    sys.stderr.write("FATAL: " + msg + "\n")
    sys.exit(1)

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
                tmin_n  = float(row[3])  # (-t) min in your CSV
                tmax_n  = float(row[4])  # (-t) max in your CSV
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
        if not isinstance(arr, list) or len(arr) != 6:
            fatal("Variable '{}' must be a list of length 6; got type {} length {}.".format(
                v, type(arr).__name__, len(arr) if isinstance(arr, list) else -1
            ))
        #endif

        for k in range(6):
            tup = arr[k]
            if (not isinstance(tup, (list, tuple))) or len(tup) != 3:
                fatal("Variable '{}' entry {} must be (tmean, value, stat).".format(v, k))
            #endif
        #endfor
    #endfor

def global_bin_to_group_and_tindex(b):
    if b < 0 or b > 23:
        fatal("Global bin out of range: {}".format(b))
    #endif
    return (b // 6), (b % 6)

def compute_migrated_values(weights, fit_map, required_varnames, renormalize, tol):
    suffixes = [
        "GEchi2FitsALUsinphi",
        "GEchi2FitsAULsinphi",
        "GEchi2FitsAULsin2phi",
        "GEchi2FitsALL",
        "GEchi2FitsALLcosphi",
    ]

    migrated = {}
    diffs = {}

    for suffix in suffixes:
        for xb_group_idx in range(4):
            xb_name = xb_group_name_from_index(xb_group_idx)
            varname = xb_name + suffix

            out_list = []
            diff_list = []

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

                delta = migrated_val - aval_target

                out_list.append([tmean_target, migrated_val, stat_target])
                diff_list.append([tmean_target, delta])
            #endfor

            migrated[varname] = out_list
            diffs[varname] = diff_list
        #endfor
    #endfor

    return migrated, diffs

def format_mathematica_list_triple(lst):
    parts = []
    for trip in lst:
        parts.append("{" + "{:.9f}".format(trip[0]) + ", " + "{:.9f}".format(trip[1]) + ", " + "{:.9f}".format(trip[2]) + "}")
    #endfor
    return "{" + ", ".join(parts) + "}"

def format_mathematica_list_pair(lst):
    parts = []
    for pair in lst:
        parts.append("{" + "{:.9f}".format(pair[0]) + ", " + "{:.9f}".format(pair[1]) + "}")
    #endfor
    return "{" + ", ".join(parts) + "}"

def write_output_files(out_path, out_diff_path, migrated, diffs, fit_map):
    suffixes = [
        "GEchi2FitsALUsinphi",
        "GEchi2FitsAULsinphi",
        "GEchi2FitsAULsin2phi",
        "GEchi2FitsALL",
        "GEchi2FitsALLcosphi",
    ]
    xb_groups = ["enpiLowxB", "enpiMidLowxB", "enpiMidHighxB", "enpiHighxB"]

    # 1) Migrated output file (triples)
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

    # 2) Diff output file:
    #    - name = {{tmean, diff}, ...};
    #    - name_ratio = {{tmean, diff/value_before}, ...};
    #
    # Fail fast if any value_before is 0 (or extremely close).
    eps = 1.0e-15

    with open(out_diff_path, "w") as f:
        for g in xb_groups:
            for s in suffixes:
                name = g + s
                if name not in diffs:
                    fatal("Internal error: diffs missing variable '{}'".format(name))
                #endif
                if name not in fit_map:
                    fatal("Internal error: fit_map missing variable '{}' needed for ratio.".format(name))
                #endif

                # Write diff line (pairs)
                rhs_diff = format_mathematica_list_pair(diffs[name])
                f.write(name + " = " + rhs_diff + ";\n")

                # Build ratio list
                ratio_list = []
                for k in range(6):
                    tmean = float(diffs[name][k][0])
                    delta = float(diffs[name][k][1])

                    before_val = float(fit_map[name][k][1])
                    if abs(before_val) < eps:
                        fatal("Cannot compute ratio for {} at index {}: before value is {:.17g} (too close to zero).".format(name, k, before_val))
                    #endif

                    ratio = delta / before_val
                    ratio_list.append([tmean, ratio])
                #endfor

                rhs_ratio = format_mathematica_list_pair(ratio_list)
                f.write(name + "_ratio = " + rhs_ratio + ";\n")
            #endfor
            f.write("\n")
        #endfor
    #endwith

def default_diff_path(out_path):
    base, ext = os.path.splitext(out_path)
    if ext == "":
        return out_path + "_diff"
    #endif
    return base + "_diff" + ext

def main():
    ap = argparse.ArgumentParser(
        description="Apply a 24x24 bin-migration matrix to selected enpi* GEchi2Fits* asymmetry lists."
    )
    ap.add_argument("migration_csv", help="CSV with 24 rows (bins 0-23) and 24 fraction columns per row.")
    ap.add_argument("fit_in_txt", help="Input Mathematica-style fit-results text file.")
    ap.add_argument("fit_out_txt", help="Output text file with migrated values (Mathematica-style).")
    ap.add_argument("--out-diff", default=None, help="Output text file for (after - before) and ratio lines. Default: <fit_out_txt with _diff inserted before extension>")
    ap.add_argument("--renormalize", action="store_true", help="Renormalize each migration row to sum to 1.0 before weighting.")
    ap.add_argument("--tol", type=float, default=1.0e-3, help="Tolerance for row-sum checks against 1.0 (strict mode).")

    args = ap.parse_args()

    out_diff = args.out_diff
    if out_diff is None:
        out_diff = default_diff_path(args.fit_out_txt)
    #endif

    weights, meta = read_migration_csv(args.migration_csv)

    fit_map = parse_fit_results_text(args.fit_in_txt)

    required = build_required_varnames()
    validate_fit_data(fit_map, required)

    migrated, diffs = compute_migrated_values(
        weights=weights,
        fit_map=fit_map,
        required_varnames=required,
        renormalize=args.renormalize,
        tol=args.tol
    )

    write_output_files(args.fit_out_txt, out_diff, migrated, diffs, fit_map)

    sys.stdout.write("Wrote migrated fit file: {}\n".format(args.fit_out_txt))
    sys.stdout.write("Wrote difference file:   {}\n".format(out_diff))

if __name__ == "__main__":
    main()
#endif