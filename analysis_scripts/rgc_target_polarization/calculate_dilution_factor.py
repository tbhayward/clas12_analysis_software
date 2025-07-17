#!/usr/bin/env python3
"""
calculate_dilution_factor.py

Module to compute dilution factor Df(xB) and its uncertainty
for each xB bin for RGC Su22 and RGC Fa22 datasets.
Saves CSV and PDF outputs, with verbose logging for progress tracking.
Includes guards against negative variance in error calculation.
"""

import numpy as np
import pandas as pd
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt
import os

# Hard-coded charge fractions (decimal) for each run period
CHARGE_FRAC = {
    "RGC_Su22": {"xNH3": 0.7239, "xC": 0.0714, "xCH2": 0.0373, "xHe": 0.0751, "xET": 0.0923},
    "RGC_Fa22": {"xNH3": 0.5849, "xC": 0.2130, "xCH2": 0.1947, "xHe": 0.0013, "xET": 0.0061},
}

def calculate_dilution_error(nNH3, nC, nCH2, nHe, nET, xNH3, xC, xCH2, xHe, xET):
    """
    Compute the uncertainty sigma_Df for the dilution factor Df,
    with guard against negative argument to sqrt.
    """
    print(f"[Error] Computing sigma_Df with counts nNH3={nNH3}, nC={nC}, nCH2={nCH2}, nHe={nHe}, nET={nET}")
    # common factor a = (-nHe*xNH3 + nNH3*xHe)
    a = -1.0 * nHe * xNH3 + nNH3 * xHe

    # compute each term
    term1 = 5438.81 * nNH3 * nET * xC**2 * xCH2**2 * xET**2 * xHe**2 \
            * a**2 \
            * (1.0 * nHe * xC * xCH2 - 1.19655 * nCH2 * xC * xHe + 0.196547 * nC * xCH2 * xHe)**2
    term2 = 85044.9 * nNH3 * nCH2 * xC**2 * xCH2**2 * xET**2 * xHe**2 \
            * a**2 \
            * (1.0 * nHe * xC * xET - 0.302592 * nET * xC * xHe - 0.697408 * nC * xET * xHe)**2
    term3 = 47470.0 * nNH3 * nC * xC**2 * xCH2**2 * xET**2 * xHe**2 \
            * a**2 \
            * (1.0 * nHe * xCH2 * xET - 0.0665285 * nET * xCH2 * xHe - 0.933472 * nCH2 * xET * xHe)**2
    term4 = 1371.83 * nHe**2 * xNH3**2 \
            * (1.0 * nHe * xC * xCH2 * xET + (-1.97748 * nET * xC * xCH2 + 6.62306 * nCH2 * xC * xET - 5.64558 * nC * xCH2 * xET) * xHe)**2 \
            * (1.0 * nHe * xC * xCH2 * xET + (0.0136533 * nET * xC * xCH2 - 1.25054 * nCH2 * xC * xET + 0.236883 * nC * xCH2 * xET) * xHe)**2
    t5_p1 = 73.2426 * nHe**2 * xNH3 * xC**2 * xCH2**2 * xET**2
    t5_p2 = nHe * xC * xCH2 * xET \
            * (2.0 * nET * xNH3 * xC * xCH2 \
               + (-183.185 * nCH2 * xNH3 * xC + 34.6998 * nC * xNH3 * xCH2 + 1.42109e-14 * nNH3 * xC * xCH2) * xET)
    t5_p3 = (-1.97748 * nET**2 * xNH3 * xC**2 * xCH2**2) \
            + nET * xC * xCH2 * (187.746 * nCH2 * xNH3 * xC - 39.9547 * nC * xNH3 * xCH2 - 145.836 * nNH3 * xC * xCH2) * xET \
            + (-606.623 * nCH2**2 * xNH3 * xC**2 \
               + nCH2 * xC * (632.002 * nC * xNH3 + 576.683 * nNH3 * xC) * xCH2 \
               + nC * (-97.9502 * nC * xNH3 - 430.847 * nNH3 * xC) * xCH2**2) * xET**2
    term5 = 0.255725 * nNH3 * nHe * (t5_p1 + t5_p2 * xHe + t5_p3 * xHe**2)**2

    denom = nNH3**3 * xHe**2 \
            * (73.2426 * nHe * xC * xCH2 * xET \
               + 1.0 * nET * xC * xCH2 * xHe \
               - 91.5925 * nCH2 * xC * xET * xHe \
               + 17.3499 * nC * xCH2 * xET * xHe)**4

    # variance must be non-negative
    var = (term1 + term2 + term3 + term4 + term5) / denom if denom != 0 else 0.0
    if var < 0:
        print(f"[Warning] Negative variance {var:.4e}; setting to zero.")
        var = 0.0

    sigma_df = 27.3473 * math.sqrt(var)
    print(f"[Error] sigma_Df = {sigma_df:.4e}")
    return sigma_df


def calculate_dilution_and_error(nNH3, nC, nCH2, nHe, nET, xNH3, xC, xCH2, xHe, xET):
    """
    Compute the dilution factor Df and its uncertainty for given counts and charge fractions.
    """
    print(f"[Dilution] Calculating Df for counts NH3={nNH3}, C={nC}, CH2={nCH2}, He={nHe}, ET={nET}")
    numerator = (27.3473 * (-1.0 * nHe * xNH3 + nNH3 * xHe)
                 * (-0.505693 * nHe * nC * nCH2 * xET
                    + (1.0 * nET * xC * xCH2 - 3.34924 * nCH2 * xC * xET + 2.85493 * nC * xCH2 * xET) * xHe))
    denominator = (nNH3 * xHe
                   * (73.2426 * nHe * xC * xCH2 * xET
                      + 1.0 * nET * xC * xCH2 * xHe
                      - 91.5925 * nCH2 * xC * xET * xHe
                      + 17.3499 * nC * xCH2 * xET * xHe))
    df = numerator / denominator if denominator != 0 else 0.0
    print(f"[Dilution] Df = {df:.4e}")
    err = calculate_dilution_error(nNH3, nC, nCH2, nHe, nET, xNH3, xC, xCH2, xHe, xET)
    print(f"[Dilution] Df error = {err:.4e}\n")
    return df, err


def _compute_period(period, trees_for_period, xB_bins):
    """
    Compute Df(xB) and sigma per xB bin for one run period.
    Returns a DataFrame with columns: x_mean, Df, Df_sigma, period.
    """
    print(f"[Period {period}] Starting dilution calculation")
    x_data = {}
    for target, tree in trees_for_period.items():
        x_data[target] = tree["x"].array(library="np")
    nbins = len(xB_bins) - 1
    rows = []
    fracs = CHARGE_FRAC[period]

    bin_idx = {t: np.digitize(x_data[t], xB_bins) - 1 for t in x_data}

    for i in range(nbins):
        mask = bin_idx["NH3"] == i
        x_mean = np.mean(x_data["NH3"][mask]) if np.any(mask) else 0.0
        nNH3 = np.sum(bin_idx["NH3"] == i)
        nC   = np.sum(bin_idx["C"]   == i)
        nCH2 = np.sum(bin_idx["CH2"] == i)
        nHe  = np.sum(bin_idx["He"]  == i)
        nET  = np.sum(bin_idx["ET"]  == i)
        print(f"[Period {period}] Bin {i+1}/{nbins}: x_mean={x_mean:.4f}, counts: NH3={nNH3}, C={nC}, CH2={nCH2}, He={nHe}, ET={nET}")
        xNH3, xC, xCH2, xHe, xET = fracs["xNH3"], fracs["xC"], fracs["xCH2"], fracs["xHe"], fracs["xET"]
        Df, Df_sigma = calculate_dilution_and_error(
            nNH3, nC, nCH2, nHe, nET,
            xNH3, xC, xCH2, xHe, xET
        )
        rows.append({"x_mean": x_mean, "Df": Df, "Df_sigma": Df_sigma})

    df = pd.DataFrame(rows)
    df["period"] = period
    print(f"[Period {period}] Completed.\n")
    return df


def calculate_and_save(trees, xB_bins):
    """
    Calculate and save dilution-factor CSV and PDF plot for Su22 and Fa22.
    """
    print("[Main] Beginning dilution-factor calculation...")
    os.makedirs("output", exist_ok=True)
    periods = ["RGC_Su22", "RGC_Fa22"]

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(_compute_period, period, trees[period], xB_bins): period for period in periods}
        results = {}
        for fut in as_completed(futures):
            period = futures[fut]
            results[period] = fut.result()
            print(f"[Main] Result received for period {period}")

    df_su = results["RGC_Su22"]
    df_fa = results["RGC_Fa22"]

    print("[Main] Computing weighted average across periods")
    w_su = 1.0 / df_su["Df_sigma"]**2
    w_fa = 1.0 / df_fa["Df_sigma"]**2
    df_avg = pd.DataFrame({
        "x_mean": df_su["x_mean"],
        "Df":     (df_su["Df"]*w_su + df_fa["Df"]*w_fa) / (w_su + w_fa),
        "Df_sigma": np.sqrt(1.0 / (w_su + w_fa)),
        "period": "weighted"
    })

    df_all = pd.concat([df_su, df_fa, df_avg], ignore_index=True)
    out_csv = "output/dilution_factor.csv"
    df_all.to_csv(out_csv, index=False)
    print(f"[Main] Saved dilution-factor CSV to {out_csv}")

    print("[Main] Generating PDF plot of Df vs xB")
    plt.errorbar(df_su["x_mean"], df_su["Df"], yerr=df_su["Df_sigma"], fmt="o", label="RGC Su22")
    plt.errorbar(df_fa["x_mean"], df_fa["Df"], yerr=df_fa["Df_sigma"], fmt="o", label="RGC Fa22")
    plt.xlabel("$x_{B}$")
    plt.ylabel("$D_{f}$")
    plt.legend(loc="upper right")
    plt.tight_layout()
    out_pdf = "output/dilution_factor.pdf"
    plt.savefig(out_pdf)
    plt.close()
    print(f"[Main] Saved dilution-factor plot to {out_pdf}\nCalculation complete.")
