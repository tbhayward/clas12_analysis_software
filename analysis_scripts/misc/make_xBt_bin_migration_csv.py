#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
make_xBt_bin_migration_csv.py

Create a 24 x 24 reconstructed-to-generated xB and -tprime bin migration CSV
for the exclusive reaction:

    e p -> e' n pi+

The output format is intentionally compatible with the older file:

    mc_bin_migration_xBt_fm.csv

Expected output format:

    #binNum,xBmin,xBmax,-tmin,-tmax,numEvents,fracfromMCbins0-23
    0,0.1,0.25,0.05,0.25,N,frac0,frac1,...,frac23
    ...
    23,0.45,0.6,1.05,1.25,N,frac0,frac1,...,frac23

The bin convention is:

    global_bin = 6 * x_bin_index + t_bin_index

where the xB bins are:

    0: 0.10 <= xB < 0.25
    1: 0.25 <= xB < 0.35
    2: 0.35 <= xB < 0.45
    3: 0.45 <= xB < 0.60

and the -tprime bins are:

    0: 0.05 <= -tprime < 0.25
    1: 0.25 <= -tprime < 0.45
    2: 0.45 <= -tprime < 0.65
    3: 0.65 <= -tprime < 0.85
    4: 0.85 <= -tprime < 1.05
    5: 1.05 <= -tprime < 1.25

For each event:

    1. Compute reconstructed t and tprime using reconstructed electron and pion.
    2. Compute generated t and tprime using generated electron and pion.
    3. Determine the reconstructed bin from reconstructed xB and reconstructed -tprime.
    4. Use the reconstructed bin to look up the Mx2 cut from the Mx2 fit parameter CSV.
    5. Require reconstructed Mx2 to satisfy:

           mu - nsigma * sigma <= Mx2 <= mu + nsigma * sigma

    6. Determine the generated bin from generated xB and generated -tprime.
    7. Fill matrix[reconstructed_bin][generated_bin].

The output row is normalized by the number of accepted events in that reconstructed bin.
Therefore each nonempty row should sum to approximately 1.

Usage example:

    python3 make_xBt_bin_migration_csv.py input.root Mx2_fit_params.csv mc_bin_migration_xBt_fm_new.csv

Example with explicit tree name:

    python3 make_xBt_bin_migration_csv.py input.root Mx2_fit_params.csv mc_bin_migration_xBt_fm_new.csv --tree PhysicsEvents

If angular branches are stored in degrees rather than radians:

    python3 make_xBt_bin_migration_csv.py input.root Mx2_fit_params.csv mc_bin_migration_xBt_fm_new.csv --angles-in-degrees

Dependencies:

    pip install uproot awkward numpy pandas
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd
import uproot


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

M_PROTON = 0.9382720813
M_NEUTRON = 0.9395654133
M_PION_PLUS = 0.13957039
M_ELECTRON = 0.00051099895

DEFAULT_BEAM_ENERGY = 10.55

XB_BINS = [
    (0.10, 0.25),
    (0.25, 0.35),
    (0.35, 0.45),
    (0.45, 0.60),
]

TPRIME_NEG_BINS = [
    (0.05, 0.25),
    (0.25, 0.45),
    (0.45, 0.65),
    (0.65, 0.85),
    (0.85, 1.05),
    (1.05, 1.25),
]

N_XB_BINS = len(XB_BINS)
N_T_BINS = len(TPRIME_NEG_BINS)
N_GLOBAL_BINS = N_XB_BINS * N_T_BINS


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

def fail(message):
    print("")
    print("ERROR:")
    print(message)
    print("")
    sys.exit(1)


def finite_mask(*arrays):
    mask = np.ones(len(arrays[0]), dtype=bool)

    for arr in arrays:
        mask &= np.isfinite(arr)
    # endfor

    return mask


def safe_sqrt(x):
    return np.sqrt(np.maximum(x, 0.0))


def lambda_kallen(a, b, c):
    return a * a + b * b + c * c - 2.0 * a * b - 2.0 * a * c - 2.0 * b * c


def find_bin(value, bins):
    for index, pair in enumerate(bins):
        low = pair[0]
        high = pair[1]

        if index == len(bins) - 1:
            if value >= low and value <= high:
                return index
            # endif
        else:
            if value >= low and value < high:
                return index
            # endif
        # endif
    # endfor

    return -1


def global_bin_index(x_index, t_index):
    if x_index < 0 or t_index < 0:
        return -1
    # endif

    if x_index >= N_XB_BINS or t_index >= N_T_BINS:
        return -1
    # endif

    return N_T_BINS * x_index + t_index


def bin_center(low, high):
    return 0.5 * (low + high)


# ---------------------------------------------------------------------------
# Four-vector and kinematic calculations
# ---------------------------------------------------------------------------

def spherical_to_cartesian(p, theta, phi, angles_in_degrees):
    if angles_in_degrees:
        theta_rad = np.deg2rad(theta)
        phi_rad = np.deg2rad(phi)
    else:
        theta_rad = theta
        phi_rad = phi
    # endif

    px = p * np.sin(theta_rad) * np.cos(phi_rad)
    py = p * np.sin(theta_rad) * np.sin(phi_rad)
    pz = p * np.cos(theta_rad)

    return px, py, pz


def minkowski_mass2(energy, px, py, pz):
    return energy * energy - px * px - py * py - pz * pz


def calculate_q_vectors(e_p, e_theta, e_phi, beam_energy, angles_in_degrees):
    e_px, e_py, e_pz = spherical_to_cartesian(
        e_p,
        e_theta,
        e_phi,
        angles_in_degrees,
    )

    e_energy = safe_sqrt(e_p * e_p + M_ELECTRON * M_ELECTRON)

    beam_energy_arr = np.full_like(e_p, beam_energy, dtype=float)
    beam_pz_arr = safe_sqrt(beam_energy_arr * beam_energy_arr - M_ELECTRON * M_ELECTRON)

    q_energy = beam_energy_arr - e_energy
    q_px = -e_px
    q_py = -e_py
    q_pz = beam_pz_arr - e_pz

    q_mass2 = minkowski_mass2(q_energy, q_px, q_py, q_pz)
    q2_positive = -q_mass2

    return q_energy, q_px, q_py, q_pz, q2_positive


def calculate_w(q_energy, q_px, q_py, q_pz):
    had_energy = M_PROTON + q_energy
    had_px = q_px
    had_py = q_py
    had_pz = q_pz

    w2 = minkowski_mass2(had_energy, had_px, had_py, had_pz)
    w = safe_sqrt(w2)

    return w, w2


def calculate_t_from_q_and_pion(
    q_energy,
    q_px,
    q_py,
    q_pz,
    pion_p,
    pion_theta,
    pion_phi,
    angles_in_degrees,
):
    pi_px, pi_py, pi_pz = spherical_to_cartesian(
        pion_p,
        pion_theta,
        pion_phi,
        angles_in_degrees,
    )

    pi_energy = safe_sqrt(pion_p * pion_p + M_PION_PLUS * M_PION_PLUS)

    diff_energy = q_energy - pi_energy
    diff_px = q_px - pi_px
    diff_py = q_py - pi_py
    diff_pz = q_pz - pi_pz

    t = minkowski_mass2(diff_energy, diff_px, diff_py, diff_pz)

    return t


def calculate_tmin_mesonic(q2_positive, w):
    """
    Compute tmin for gamma* p -> pi+ n.

    This uses the mesonic definition:

        t = (q - p_pi)^2

    In the gamma* p center-of-mass frame:

        t = -Q2 + m_pi^2 - 2 * (E_gamma_star * E_pi_star
            - |q_star| * |p_pi_star| * cos(theta_star))

    The minimum |t|, usually called tmin, corresponds to forward meson
    production, cos(theta_star) = +1.

    For physical events this value is normally negative and closest to zero.
    """

    w2 = w * w

    e_gamma_star = (w2 - M_PROTON * M_PROTON - q2_positive) / (2.0 * w)
    q_star = safe_sqrt(e_gamma_star * e_gamma_star + q2_positive)

    e_pi_star = (w2 + M_PION_PLUS * M_PION_PLUS - M_NEUTRON * M_NEUTRON) / (2.0 * w)

    lam = lambda_kallen(w2, M_PION_PLUS * M_PION_PLUS, M_NEUTRON * M_NEUTRON)
    p_pi_star = safe_sqrt(lam) / (2.0 * w)

    tmin = (
        -q2_positive
        + M_PION_PLUS * M_PION_PLUS
        - 2.0 * (e_gamma_star * e_pi_star - q_star * p_pi_star)
    )

    return tmin


def calculate_tprime_for_sample(
    e_p,
    e_theta,
    e_phi,
    pion_p,
    pion_theta,
    pion_phi,
    beam_energy,
    angles_in_degrees,
):
    q_energy, q_px, q_py, q_pz, q2_positive = calculate_q_vectors(
        e_p,
        e_theta,
        e_phi,
        beam_energy,
        angles_in_degrees,
    )

    w, w2 = calculate_w(q_energy, q_px, q_py, q_pz)

    t = calculate_t_from_q_and_pion(
        q_energy,
        q_px,
        q_py,
        q_pz,
        pion_p,
        pion_theta,
        pion_phi,
        angles_in_degrees,
    )

    tmin = calculate_tmin_mesonic(q2_positive, w)

    tprime = t - tmin

    return t, tmin, tprime, q2_positive, w, w2


# ---------------------------------------------------------------------------
# Mx2 cut table
# ---------------------------------------------------------------------------

def load_mx2_cut_table(mx2_csv_path):
    required_columns = [
        "slice_name",
        "x_min",
        "x_max",
        "t_min",
        "t_max",
        "mu",
        "sigma",
        "n_entries",
        "fit_success",
    ]

    if not os.path.exists(mx2_csv_path):
        fail("Mx2 fit parameter CSV does not exist: {}".format(mx2_csv_path))
    # endif

    df = pd.read_csv(mx2_csv_path)

    missing_columns = []
    for col in required_columns:
        if col not in df.columns:
            missing_columns.append(col)
        # endif
    # endfor

    if len(missing_columns) > 0:
        fail(
            "Mx2 fit parameter CSV is missing required columns:\n{}".format(
                "\n".join(missing_columns)
            )
        )
    # endif

    if len(df) != N_GLOBAL_BINS:
        fail(
            "Expected {} Mx2 cut rows, but found {} rows in {}".format(
                N_GLOBAL_BINS,
                len(df),
                mx2_csv_path,
            )
        )
    # endif

    cuts = {}

    for row_index, row in df.iterrows():
        x_min = float(row["x_min"])
        x_max = float(row["x_max"])
        t_min = float(row["t_min"])
        t_max = float(row["t_max"])
        mu = float(row["mu"])
        sigma = float(row["sigma"])
        fit_success = int(row["fit_success"])

        if fit_success != 1:
            fail(
                "Mx2 fit row {} has fit_success != 1. Refusing to use it.".format(
                    row_index
                )
            )
        # endif

        x_index = find_bin(bin_center(x_min, x_max), XB_BINS)
        t_index = find_bin(bin_center(t_min, t_max), TPRIME_NEG_BINS)
        gbin = global_bin_index(x_index, t_index)

        if gbin < 0:
            fail(
                "Could not map Mx2 cut row {} to an analysis bin: "
                "x=({}, {}), -tprime=({}, {})".format(
                    row_index,
                    x_min,
                    x_max,
                    t_min,
                    t_max,
                )
            )
        # endif

        if gbin in cuts:
            fail("Duplicate Mx2 cut entry for global bin {}".format(gbin))
        # endif

        cuts[gbin] = {
            "mu": mu,
            "sigma": sigma,
            "x_min": x_min,
            "x_max": x_max,
            "t_min": t_min,
            "t_max": t_max,
        }
    # endfor

    missing_bins = []
    for gbin in range(N_GLOBAL_BINS):
        if gbin not in cuts:
            missing_bins.append(gbin)
        # endif
    # endfor

    if len(missing_bins) > 0:
        fail(
            "Mx2 cut table is missing these global bins:\n{}".format(
                ", ".join(str(x) for x in missing_bins)
            )
        )
    # endif

    return cuts


# ---------------------------------------------------------------------------
# ROOT reading
# ---------------------------------------------------------------------------

def required_branches():
    branches = [
        "e_p",
        "mc_e_p",
        "e_theta",
        "mc_e_theta",
        "e_phi",
        "mc_e_phi",
        "p_p",
        "mc_p_p",
        "p_theta",
        "mc_p_theta",
        "p_phi",
        "mc_p_phi",
        "x",
        "mc_x",
        "Mx2",
    ]

    return branches


def open_tree(input_root, tree_name):
    if not os.path.exists(input_root):
        fail("Input ROOT file does not exist: {}".format(input_root))
    # endif

    try:
        root_file = uproot.open(input_root)
    except Exception as exc:
        fail("Could not open ROOT file {}\n{}".format(input_root, exc))
    # endtry

    if tree_name not in root_file:
        available_keys = list(root_file.keys())
        fail(
            "Tree '{}' not found in {}\nAvailable keys:\n{}".format(
                tree_name,
                input_root,
                "\n".join(str(k) for k in available_keys),
            )
        )
    # endif

    return root_file[tree_name]


def validate_branches(tree, branches):
    available = set(tree.keys())

    missing = []
    for branch in branches:
        if branch not in available:
            missing.append(branch)
        # endif
    # endfor

    if len(missing) > 0:
        fail(
            "Input tree is missing required branches:\n{}".format(
                "\n".join(missing)
            )
        )
    # endif


# ---------------------------------------------------------------------------
# Main matrix construction
# ---------------------------------------------------------------------------

def build_migration_matrix(
    input_root,
    tree_name,
    mx2_csv_path,
    output_csv_path,
    beam_energy,
    angles_in_degrees,
    nsigma,
    batch_size,
):
    mx2_cuts = load_mx2_cut_table(mx2_csv_path)

    tree = open_tree(input_root, tree_name)
    branches = required_branches()
    validate_branches(tree, branches)

    matrix_counts = np.zeros((N_GLOBAL_BINS, N_GLOBAL_BINS), dtype=np.float64)
    row_denominators = np.zeros(N_GLOBAL_BINS, dtype=np.float64)

    total_events = 0
    total_finite_events = 0
    total_reco_in_range = 0
    total_pass_mx2 = 0
    total_mc_in_range = 0
    total_filled = 0

    print("")
    print("Building xB and -tprime migration matrix")
    print("  input ROOT file  : {}".format(input_root))
    print("  tree             : {}".format(tree_name))
    print("  Mx2 cut CSV      : {}".format(mx2_csv_path))
    print("  output CSV       : {}".format(output_csv_path))
    print("  beam energy      : {:.6f} GeV".format(beam_energy))
    print("  pion branches    : p_p, p_theta, p_phi")
    print("  angles in degrees: {}".format(angles_in_degrees))
    print("  Mx2 cut width    : mu +/- {:.3f} sigma".format(nsigma))
    print("")

    for arrays in tree.iterate(
        branches,
        step_size=batch_size,
        library="np",
    ):
        n_batch = len(arrays["x"])
        total_events += n_batch

        reco_t, reco_tmin, reco_tprime, reco_q2_calc, reco_w_calc, reco_w2_calc = calculate_tprime_for_sample(
            arrays["e_p"].astype(float),
            arrays["e_theta"].astype(float),
            arrays["e_phi"].astype(float),
            arrays["p_p"].astype(float),
            arrays["p_theta"].astype(float),
            arrays["p_phi"].astype(float),
            beam_energy,
            angles_in_degrees,
        )

        mc_t, mc_tmin, mc_tprime, mc_q2_calc, mc_w_calc, mc_w2_calc = calculate_tprime_for_sample(
            arrays["mc_e_p"].astype(float),
            arrays["mc_e_theta"].astype(float),
            arrays["mc_e_phi"].astype(float),
            arrays["mc_p_p"].astype(float),
            arrays["mc_p_theta"].astype(float),
            arrays["mc_p_phi"].astype(float),
            beam_energy,
            angles_in_degrees,
        )

        reco_x = arrays["x"].astype(float)
        mc_x = arrays["mc_x"].astype(float)
        reco_mx2 = arrays["Mx2"].astype(float)

        reco_neg_tprime = -reco_tprime
        mc_neg_tprime = -mc_tprime

        mask_finite = finite_mask(
            reco_x,
            mc_x,
            reco_neg_tprime,
            mc_neg_tprime,
            reco_mx2,
            reco_t,
            reco_tmin,
            mc_t,
            mc_tmin,
        )

        total_finite_events += int(np.count_nonzero(mask_finite))

        finite_indices = np.where(mask_finite)[0]

        for i in finite_indices:
            reco_x_index = find_bin(reco_x[i], XB_BINS)
            reco_t_index = find_bin(reco_neg_tprime[i], TPRIME_NEG_BINS)
            reco_bin = global_bin_index(reco_x_index, reco_t_index)

            if reco_bin < 0:
                continue
            # endif

            total_reco_in_range += 1

            mx2_cut = mx2_cuts[reco_bin]
            mx2_low = mx2_cut["mu"] - nsigma * mx2_cut["sigma"]
            mx2_high = mx2_cut["mu"] + nsigma * mx2_cut["sigma"]

            if reco_mx2[i] < mx2_low or reco_mx2[i] > mx2_high:
                continue
            # endif

            total_pass_mx2 += 1

            mc_x_index = find_bin(mc_x[i], XB_BINS)
            mc_t_index = find_bin(mc_neg_tprime[i], TPRIME_NEG_BINS)
            mc_bin = global_bin_index(mc_x_index, mc_t_index)

            if mc_bin < 0:
                continue
            # endif

            total_mc_in_range += 1

            matrix_counts[reco_bin, mc_bin] += 1.0
            row_denominators[reco_bin] += 1.0
            total_filled += 1
        # endfor
    # endfor

    write_migration_csv(
        output_csv_path,
        matrix_counts,
        row_denominators,
    )

    print("")
    print("Done.")
    print("Summary:")
    print("  total events read                         : {}".format(total_events))
    print("  finite events                             : {}".format(total_finite_events))
    print("  reco events inside xB/-tprime bins        : {}".format(total_reco_in_range))
    print("  reco events passing bin-dependent Mx2 cut : {}".format(total_pass_mx2))
    print("  events also inside generated bins         : {}".format(total_mc_in_range))
    print("  events filled into migration matrix       : {}".format(total_filled))
    print("")
    print("Per-row event counts:")
    for gbin in range(N_GLOBAL_BINS):
        x_index = gbin // N_T_BINS
        t_index = gbin % N_T_BINS

        x_low, x_high = XB_BINS[x_index]
        t_low, t_high = TPRIME_NEG_BINS[t_index]

        print(
            "  bin {:2d}: xB [{:.2f}, {:.2f}], -tprime [{:.2f}, {:.2f}] -> {}".format(
                gbin,
                x_low,
                x_high,
                t_low,
                t_high,
                int(row_denominators[gbin]),
            )
        )
    # endfor

    print("")
    print("Wrote:")
    print("  {}".format(output_csv_path))
    print("")


def write_migration_csv(output_csv_path, matrix_counts, row_denominators):
    output_dir = os.path.dirname(os.path.abspath(output_csv_path))

    if output_dir != "" and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # endif

    with open(output_csv_path, "w", newline="") as outfile:
        writer = csv.writer(outfile)

        outfile.write("#binNum,xBmin,xBmax,-tmin,-tmax,numEvents,fracfromMCbins0-23\n")

        for reco_bin in range(N_GLOBAL_BINS):
            x_index = reco_bin // N_T_BINS
            t_index = reco_bin % N_T_BINS

            x_low, x_high = XB_BINS[x_index]
            t_low, t_high = TPRIME_NEG_BINS[t_index]

            denominator = row_denominators[reco_bin]

            if denominator > 0.0:
                fractions = matrix_counts[reco_bin, :] / denominator
            else:
                fractions = np.zeros(N_GLOBAL_BINS, dtype=float)
            # endif

            row = [
                reco_bin,
                "{:.10g}".format(x_low),
                "{:.10g}".format(x_high),
                "{:.10g}".format(t_low),
                "{:.10g}".format(t_high),
                int(denominator),
            ]

            for frac in fractions:
                row.append("{:.10g}".format(float(frac)))
            # endfor

            writer.writerow(row)
        # endfor
    # endwith


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a 24 x 24 xB/-tprime migration CSV from a ROOT tree "
            "containing reconstructed and generated kinematics."
        )
    )

    parser.add_argument(
        "input_root",
        help="Input ROOT file containing the tree.",
    )

    parser.add_argument(
        "mx2_fit_params_csv",
        help="CSV containing bin-dependent Mx2 fit parameters.",
    )

    parser.add_argument(
        "output_csv",
        help="Output migration CSV path.",
    )

    parser.add_argument(
        "--tree",
        default="PhysicsEvents",
        help="Tree name. Default: PhysicsEvents",
    )

    parser.add_argument(
        "--beam-energy",
        type=float,
        default=DEFAULT_BEAM_ENERGY,
        help="Input lepton beam energy in GeV. Default: 10.55",
    )

    parser.add_argument(
        "--angles-in-degrees",
        action="store_true",
        help=(
            "Use this if e_theta, e_phi, p_theta, p_phi, etc. are stored "
            "in degrees. By default the script assumes radians."
        ),
    )

    parser.add_argument(
        "--nsigma",
        type=float,
        default=3.0,
        help="Mx2 cut width: mu +/- nsigma*sigma. Default: 3.0",
    )

    parser.add_argument(
        "--batch-size",
        default="200 MB",
        help="uproot iterate step size. Default: 200 MB",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    build_migration_matrix(
        input_root=args.input_root,
        tree_name=args.tree,
        mx2_csv_path=args.mx2_fit_params_csv,
        output_csv_path=args.output_csv,
        beam_energy=args.beam_energy,
        angles_in_degrees=args.angles_in_degrees,
        nsigma=args.nsigma,
        batch_size=args.batch_size,
    )


if __name__ == "__main__":
    main()
# endif